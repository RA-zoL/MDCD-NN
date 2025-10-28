#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, socket, pickle
import numpy as NP
import ase
from ase.data import atomic_numbers
from ase.calculators.calculator import Calculator, all_changes

try:
    import torch
    from deepmd.infer import DeepPot as DPEval
    from deepmd.calculator import DP as DPCalculator
    HAS_DPMD = True
except ImportError:
    HAS_DPMD = False

class DPEnsemble(Calculator):
    '''
    ASE calculator for DeePMD NNP ensemble
    '''
    
    name = "DPEnsemble"
    implemented_properties = ["energy", "forces", "dipole"] # model deviation -> dipole

    @staticmethod
    def default_fparam(atoms):
        if "charge" not in atoms.info:
            atoms.info["charge"] = 0
        if "spin" not in atoms.info:
            atoms.info["spin"] = 1
        return NP.array([[atoms.info["charge"], atoms.info["spin"]]])

    def __init__(self, Models = None, TypeMap = None, get_fparam = None, **kwargs):
        '''
        Models: list[str], DeePMD model paths
        TypeMap: dict[int, int], type map from atomic number to DeePMD atomic type
        get_fparam: callable, to get fparam for atoms
        '''
        super().__init__(**kwargs)
        self.Models = []
        for M in Models:
            try:
                self.Models.append(DPEval(M))
            except:
                print(f'Warning: DPMD model {M} cannot be loaded!')
        if len(self.Models) == 0:
            print(f'Error: No models are loaded.')
            return
        self.TypeMap = dict(TypeMap)
        self.get_fparam = get_fparam if get_fparam is not None else self.default_fparam

    def calculate(self, atoms = None, properties = ["energy", "forces"], system_changes = all_changes):
        '''
        results:
            energy, float in eV, forces, float(natoms * 3) in eV/A
            E_devi, float in eV/atom, F_devi, float in eV/A
        '''
        super().calculate(atoms, properties, system_changes)
        Type = NP.array([self.TypeMap[i] for i in atoms.numbers])
        Coord = atoms.positions.reshape((1, len(atoms), 3))
        Cell = None
        FParam = self.get_fparam(atoms)
        Es = []
        Fs = []
        for i in range(len(self.Models)):
            CurE, CurF, CurV = self.Models[i].eval(Coord, Cell, Type, atomic = False, fparam = FParam)
            Es.append(CurE)
            Fs.append(CurF)
        Es = NP.array(Es).reshape([len(self.Models)])
        Fs = NP.array(Fs).reshape([len(self.Models), len(atoms), 3])
        PredE = NP.average(Es)
        DeviE = NP.std(Es) / len(atoms)
        PredF = NP.average(Fs, axis = 0)
        DeviF = NP.max(NP.linalg.norm(NP.std(Fs, axis = 0), axis = -1))
        self.results = {"energy": PredE, "forces": PredF, "E_devi": DeviE, "F_devi": DeviF}
        atoms.info["F_devi"] = DeviF # save model deviation in atoms

def calc_hessian(A: ase.Atoms, HStep = 0.001):
    '''
    hessian calculator for inner invocation
    HStep: float, geometry displacement for numerical hessian
    '''
    NAtoms = len(A)
    H = NP.zeros((NAtoms * 3, NAtoms * 3), dtype = 'd')
    OriginalPos = A.positions
    LeftC = NP.array([OriginalPos for i in range(NAtoms * 3)])
    RightC = NP.array([OriginalPos for i in range(NAtoms * 3)])
    for i in range(NAtoms):
        for j in range(3):
            LeftC[i * 3 + j][i][j] -= HStep
            RightC[i * 3 + j][i][j] += HStep
    LeftG = NP.zeros((NAtoms * 3, NAtoms * 3), dtype = 'd')
    RightG = NP.zeros((NAtoms * 3, NAtoms * 3), dtype = 'd')
    for i in range(NAtoms * 3):
        A.positions = LeftC[i]
        A.calc.calculate(atoms = A, properties = ["energy", "forces"], system_changes = all_changes)
        LeftG[i] = -A.calc.results["forces"].reshape([-1])
        A.positions = RightC[i]
        A.calc.calculate(atoms = A, properties = ["energy", "forces"], system_changes = all_changes)
        RightG[i] = -A.calc.results["forces"].reshape([-1])
    for i in range(NAtoms * 3):
        for j in range(i, NAtoms * 3):
            H[j][i] = H[i][j] = (RightG[i][j] + RightG[j][i] - LeftG[i][j] - LeftG[j][i]) / 4.0 / HStep
    A.positions = OriginalPos
    return H

def ServerMain():
    '''
    Predicting energies and forces externally
    Usage: DPA3_Server.py model_1 model_2 ...
    '''
    
    # Server initialization
    MaxBufferSize = 131072 # should be able to receive >500 atoms once
    TerminatingCode = -1
    ErrorCode = 5
    Host = socket.gethostname()
    Port = 31079 # A prime number brings good luck
    ServerSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    ServerSocket.bind((Host, Port))
    ServerSocket.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 128 * MaxBufferSize)
    ServerSocket.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, 128 * MaxBufferSize)
    Termination = False
    
    # DPA-3 model initialization
    Models = sys.argv[1:]
    TypeMap = {1:0, 6:1, 7:2, 8:3} # H C N O
    Calc = DPEnsemble(Models = Models, TypeMap = TypeMap)
    HessianStep = 0.001 # A
    
    MinDeviF = 0.1 # eV/A
    MaxDeviF = 0.5 # eV/A
    Acc = Candi = Fail = 0
    CountE = CountH = 0
    print(f'DPA-3 server initialized on {Host}:{Port}')

    while not Termination:
        # receiving data
        Data, ClientAddress = ServerSocket.recvfrom(MaxBufferSize)
        try:
            Data = pickle.loads(Data)
            if not Data:
                ServerSocket.sendto(pickle.dumps(ErrorCode), ClientAddress)
            if Data == TerminatingCode:
                Termination = True
                ServerSocket.sendto(pickle.dumps(TerminatingCode), ClientAddress)
                break

            # calculate
            FilePrefix, Derives, Charge, Spin = Data
            CurNumbers = NP.load(f'{FilePrefix}.atomic_number.npy')
            CurCoord = NP.load(f'{FilePrefix}.coord.npy')
            CurCfm = ase.Atoms(numbers = CurNumbers, positions = CurCoord)
            CurCfm.info["charge"] = Charge
            CurCfm.info["spin"] = Spin
            CurCfm.calc = Calc
            CurCfm.calc.calculate(atoms = CurCfm)

            PredE = Calc.results["energy"]
            DeviE = Calc.results["E_devi"]
            NP.save(f'''{FilePrefix}.force.npy''', Calc.results["forces"])
            DeviF = Calc.results["F_devi"]
            CountE += 1

            # Numerical hessian
            if Derives > 1:
                PredH = calc_hessian(CurCfm, HessianStep)
                NP.save(f'''{FilePrefix}.hessian.npy''', PredH)
                CountH += 1
            else:
                PredH = 0
            
            # sending results
            ServerSocket.sendto(pickle.dumps((PredE, DeviE, DeviF)), ClientAddress)

            # statistics
            if DeviF < MinDeviF:
                Acc += 1
            elif DeviF < MaxDeviF:
                Candi += 1
            else:
                Fail += 1
        except:
            ServerSocket.sendto(pickle.dumps(ErrorCode), ClientAddress)

    ServerSocket.close()
    
    print(f'''ML predicted engrads: {CountE}''')
    print(f'''ML predicted hessians: {CountH}''')
    print(f'''ML force deviation threshold: Min= {MinDeviF:.2f} Max= {MaxDeviF:.2f}''')
    print(f'''ML accurate predictions: {Acc}''')
    print(f'''ML candidate predictions: {Candi}''')
    print(f'''ML failed predictions: {Fail}''')

if __name__ == "__main__":
    ServerMain()

