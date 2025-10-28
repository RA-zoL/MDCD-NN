#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, socket, pickle
import numpy as NP
import ase.units

MAXBUFFERSIZE = 131072 # should be able to receive >500 atoms
ERROR_CODE = 5

# This is an interface invoked by Gaussian 16 external 

# Totally 6 subsections
# 0. process command line
# write the route line in keyword like: external='DPA3_Client.py'
# command line invoked by Gaussian should be DPA3_Client.py R (Gaussian external args)
ServerHost = socket.gethostname()
ServerPort = 31079 # A prime number brings good luck
SplitTag = sys.argv.index("R")
InputFileName, OutputFileName, MsgFileName = sys.argv[SplitTag + 1: SplitTag + 4]
ScratchDir = InputFileName[:InputFileName.find("Gau-") - 1]
PID = InputFileName[InputFileName.find("Gau-") + 4:InputFileName.find(".E")]

# 1. Client initialization
ClientSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
ClientSocket.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, MAXBUFFERSIZE)
ClientSocket.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, MAXBUFFERSIZE)

# 2. Read Gaussian external inputs
with open(InputFileName, 'r') as InputFile:
    Text = InputFile.readlines()
splitline = (' '.join(Text[0].strip('\n').split())).split()
NAtoms = int(splitline[0])
Derives = int(splitline[1])
Charge = int(splitline[2])
Spin = int(splitline[3])

AtomicNumbers = NP.zeros(NAtoms, dtype = 'i')
Positions = NP.zeros((NAtoms, 3), dtype = 'd')
for i in range(NAtoms):
    splitline = (' '.join(Text[i + 1].strip('\n').split())).split()
    AtomicNumbers[i] = int(splitline[0])
    Positions[i][0] = float(splitline[1])
    Positions[i][1] = float(splitline[2])
    Positions[i][2] = float(splitline[3])
Positions *= ase.units.Bohr

# 3. Send atomic numbers, positions, total charge and total spin to server
FilePrefix = f'''{ScratchDir}/{PID}'''
NP.save(f'''{FilePrefix}.atomic_number.npy''', AtomicNumbers)
NP.save(f'''{FilePrefix}.coord.npy''', Positions)
Sent = ClientSocket.sendto(pickle.dumps((FilePrefix, Derives, Charge, Spin)), (ServerHost, ServerPort))

# 4. Receive results from server
Data, ServerAddress = ClientSocket.recvfrom(MAXBUFFERSIZE)
try:
    Data = pickle.loads(Data)
    if not Data:
        exit(1)
    if Data == ERROR_CODE:
        exit(1)
    Energy = Data[0] / ase.units.Hartree
    Gradients = -NP.load(f'''{FilePrefix}.force.npy''').reshape([NAtoms, 3]) / ase.units.Hartree * ase.units.Bohr
    ierr = os.system(f'''rm {FilePrefix}.atomic_number.npy {FilePrefix}.coord.npy {FilePrefix}.force.npy''')
    if Derives > 1:
        Hessians = NP.load(f'''{FilePrefix}.hessian.npy''') / ase.units.Hartree * ase.units.Bohr ** 2
        ierr = os.system(f'''rm {FilePrefix}.hessian.npy''')
    DeviE = Data[1]
    DeviF = Data[2]
except:
    ClientSocket.close()

# 5. Write Gaussian external outputs
with open(MsgFileName, 'w') as Log:
    Log.write(f'''SCF Done:  E(DPA-3) = {Energy:16.9f}     A.U. after 0 cycles\n''')
    Log.write(f'''Deviation: devi_E = {DeviE * 1.0e3:8.4f} meV/NAtoms  devi_F = {DeviF * 1.0e3:8.4f} meV/A\n''')

with open(OutputFileName, 'w') as OutputFile:
    # energy
    OutputFile.write(f'''{Energy:20.12f}{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n''')
    if Derives > 0:
        # gradients
        for i in range(NAtoms):
            OutputFile.write(f'''{Gradients[i][0]:20.12f}'''
                f'''{Gradients[i][1]:20.12f}{Gradients[i][2]:20.12f}\n''')
    if Derives > 1:
        # polarizability
        OutputFile.write(f'''{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n''')
        OutputFile.write(f'''{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n''')
        # dipole derivatives
        for i in range(3 * NAtoms):
            OutputFile.write(f'''{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n''')
        # hessian
        cnt = 0
        for i in range(3 * NAtoms):
            for j in range(i + 1):
                OutputFile.write(f'''{Hessians[i][j]:20.12f}''')
                cnt += 1
                if cnt == 3:
                    cnt = 0
                    OutputFile.write(f'''\n''')
        OutputFile.write(f'''\n''')

# 6. Finalization
ClientSocket.close()
