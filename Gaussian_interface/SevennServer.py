#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, socket, pickle
import numpy as NP
import ase, ase.units, ase.io, ase.calculators.calculator
import torch
from sevenn.sevennet_calculator import SevenNetCalculator

def SevennServerMain():
    '''
    Predicting energies and forces externally
    '''
    
    # Server initialization
    MaxBufferSize = 131072 # should be able to receive >500 atoms once
    TerminatingCode = -1
    ErrorCode = -5
    Host = socket.gethostname()
    Port = 31079 # A prime number brings good luck
    ServerSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    ServerSocket.bind((Host, Port))
    ServerSocket.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 128 * MaxBufferSize)
    ServerSocket.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, 128 * MaxBufferSize)
    Terminate = False
    
    # Sevenn model initialization
    NModels = int(sys.argv[1])
    M = [SevenNetCalculator(sys.argv[i + 2], device = "cuda" if torch.cuda.is_available() else "cpu") for i in range(NModels)]
    HessianStep = 0.001 # A
    CurE = [None for i in range(NModels)]
    CurF = [None for i in range(NModels)]
    MinDeviF = 0.1 # eV/A
    MaxDeviF = 0.6 # eV/A
    Acc = Candi = Fail = 0
    CountE = CountH = 0
    ierr = os.system(f'''touch {sys.argv[NModels + 2]}.initialized''')

    while not Terminate:
        # receiving data
        Data, ClientAddress = ServerSocket.recvfrom(MaxBufferSize)
        if True:
            Data = pickle.loads(Data)
            if not Data:
                ServerSocket.sendto(pickle.dumps(ErrorCode), ClientAddress)
            if Data == TerminatingCode:
                Terminate = True
                ServerSocket.sendto(pickle.dumps(TerminatingCode), ClientAddress)
                break

            # Applying Sevenn ensemble
            FilePrefix = Data[0]
            Derives = Data[1]
            CurType = NP.load(f'''{FilePrefix}.atomic_number.npy''')
            NAtoms = len(CurType)
            CurCoord = NP.load(f'''{FilePrefix}.coord.npy''')
            A = ase.Atoms(numbers = CurType, positions = CurCoord)

            for i in range(NModels):
                A.calc = M[i]
                M[i].calculate(atoms = A, properties = ["energy", "forces"],
                    system_changes = ase.calculators.calculator.all_changes)
                CurE[i] = M[i].results["energy"]
                CurF[i] = M[i].results["forces"]
            CountE += 1

            PredE = NP.average(NP.array(CurE))
            DeviE = (NP.sum((NP.array(CurE) - PredE) ** 2) / NModels) ** 0.5 / NAtoms
            PredFs = NP.array(CurF)
            PredF = NP.average(PredFs, axis = 0)
            NP.save(f'''{FilePrefix}.force.npy''', PredF)
            DeviF = NP.max(NP.linalg.norm(NP.std(PredFs, axis = 0), axis = -1))

            # Numerical hessian
            if Derives > 1:
                PredH = NP.zeros(NAtoms * NAtoms * 9, dtype = 'd').reshape([NAtoms * 3, NAtoms * 3])
                LeftC = NP.array([CurCoord for i in range(NAtoms * 3)])
                RightC = NP.array([CurCoord for i in range(NAtoms * 3)])
                for i in range(NAtoms):
                    for j in range(3):
                        LeftC[i * 3 + j][i][j] -= HessianStep
                        RightC[i * 3 + j][i][j] += HessianStep
                LeftG = NP.zeros(NAtoms * 3 * NAtoms * 3, dtype = 'd').reshape([NAtoms * 3, NAtoms * 3])
                RightG = NP.zeros(NAtoms * 3 * NAtoms * 3, dtype = 'd').reshape([NAtoms * 3, NAtoms * 3])
                for i in range(NAtoms):
                    for j in range(3):
                        for k in range(NModels):
                            B = ase.Atoms(numbers = CurType, positions = LeftC[i * 3 + j])
                            B.calc = M[k]
                            M[k].calculate(atoms = B, properties = ["energy", "forces"],
                                system_changes = ase.calculators.calculator.all_changes)
                            CurE[k] = M[k].results["energy"]
                            CurF[k] = M[k].results["forces"]
                        LeftG[i * 3 + j] = -NP.average(NP.array(CurF), axis = 0).reshape([-1])
                        for k in range(NModels):
                            B = ase.Atoms(numbers = CurType, positions = RightC[i * 3 + j])
                            B.calc = M[k]
                            M[k].calculate(atoms = B, properties = ["energy", "forces"],
                                system_changes = ase.calculators.calculator.all_changes)
                            CurE[k] = M[k].results["energy"]
                            CurF[k] = M[k].results["forces"]
                        RightG[i * 3 + j] = -NP.average(NP.array(CurF), axis = 0).reshape([-1])
                for i in range(NAtoms * 3):
                    for j in range(i, NAtoms * 3):
                        PredH[j][i] = PredH[i][j] = (RightG[i][j] + RightG[j][i] - LeftG[i][j] - LeftG[j][i]) / 4 / HessianStep
                CountH += 1
                NP.save(f'''{FilePrefix}.hessian.npy''', PredH)
            else:
                PredH = 0
            ServerSocket.sendto(pickle.dumps((PredE, DeviE, DeviF)), ClientAddress)
            if DeviF < MinDeviF:
                Acc += 1
            elif DeviF < MaxDeviF:
                Candi += 1
            else:
                Fail += 1
        else:
            ServerSocket.sendto(pickle.dumps(ErrorCode), ClientAddress)

    ServerSocket.close()
    
    print(f'''ML predicted engrads: {CountE}''')
    print(f'''ML predicted hessians: {CountH}''')
    print(f'''ML force deviation threshold: Min= {MinDeviF:.2f} Max= {MaxDeviF:.2f}''')
    print(f'''ML accurate predictions: {Acc}''')
    print(f'''ML candidate predictions: {Candi}''')
    print(f'''ML failed predictions: {Fail}''')

if __name__ == "__main__":
    SevennServerMain()

