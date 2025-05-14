#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, socket, pickle
import numpy as NP
import ase.units

# process command line
SplitTag = 1
while sys.argv[SplitTag] != "R":
	SplitTag += 1
InputFileName, OutputFileName, MsgFileName = sys.argv[SplitTag + 1: SplitTag + 4]
ScratchDir = InputFileName[:InputFileName.find("Gau-") - 1]
PID = InputFileName[InputFileName.find("Gau-") + 4:InputFileName.find(".E")]

# Server initialization
MaxBufferSize = 16384
ErrorCode = -5
Host = socket.gethostname()
Port = 31079 # A prime number brings good luck
ClientSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
ClientSocket.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 2 * MaxBufferSize)
ClientSocket.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, 2 * MaxBufferSize)

# Reading Gaussian external inputs
with open(InputFileName, 'r') as InputFile:
    Text = InputFile.readlines()
splitline = (' '.join(Text[0].strip('\n').split())).split()
NAtoms = int(splitline[0])
Derives = int(splitline[1])
Charge = int(splitline[2])
Spin = int(splitline[3])

AtomicNumbers = NP.zeros(NAtoms, dtype = 'i')
Positions = NP.zeros(3 * NAtoms, dtype = 'd').reshape([NAtoms, 3])
for i in range(NAtoms):
    splitline = (' '.join(Text[i + 1].strip('\n').split())).split()
    AtomicNumbers[i] = int(splitline[0])
    Positions[i][0] = float(splitline[1])
    Positions[i][1] = float(splitline[2])
    Positions[i][2] = float(splitline[3])
Positions *= ase.units.Bohr

# Predicting energies and forces externally
FilePrefix = f'''{ScratchDir}/{PID}'''
NP.save(f'''{FilePrefix}.atomic_number.npy''', AtomicNumbers)
NP.save(f'''{FilePrefix}.coord.npy''', Positions)
Sent = ClientSocket.sendto(pickle.dumps((FilePrefix, Derives)), (Host, Port))
Data, ServerAddress = ClientSocket.recvfrom(MaxBufferSize)
if True:
    Data = pickle.loads(Data)
    if not Data:
        exit(1)
    if Data == ErrorCode:
        exit(1)
    Energy = Data[0] / ase.units.Hartree
    Gradients = -NP.load(f'''{FilePrefix}.force.npy''') / ase.units.Hartree * ase.units.Bohr
    ierr = os.system(f'''rm {FilePrefix}.atomic_number.npy {FilePrefix}.coord.npy {FilePrefix}.force.npy''')
    if Derives > 1:
        Hessians = NP.load(f'''{FilePrefix}.hessian.npy''') / ase.units.Hartree * ase.units.Bohr ** 2
        ierr = os.system(f'''rm {FilePrefix}.hessian.npy''')
    DeviE = Data[1]
    DeviF = Data[2]
else:
    ClientSocket.close()

with open(MsgFileName, 'w') as Log:
    Log.write(f'''SCF Done:  E(SeveNN) = {Energy:16.9f}     A.U. after 0 cycles\n''')
    Log.write(f'''Deviation: devi_E = {DeviE * 1000.0:8.4f} meV/NAtoms  devi_F = {DeviF * 1000.0:8.4f} meV/A\n''')

with open(OutputFileName, 'w') as OutputFile:
    OutputFile.write(f'''{Energy:20.12f}{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n''')
    if Derives > 0:
        for i in range(NAtoms):
            OutputFile.write(f'''{Gradients[i][0]:20.12f}'''
                f'''{Gradients[i][1]:20.12f}{Gradients[i][2]:20.12f}\n''')
    if Derives > 1:
        OutputFile.write(f'''{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n''')
        OutputFile.write(f'''{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n''')
        for i in range(3 * NAtoms):
            OutputFile.write(f'''{0.0:20.12f}{0.0:20.12f}{0.0:20.12f}\n''')
        cnt = 0
        for i in range(3 * NAtoms):
            for j in range(i + 1):
                OutputFile.write(f'''{Hessians[i][j]:20.12f}''')
                cnt += 1
                if cnt == 3:
                    cnt = 0
                    OutputFile.write(f'''\n''')
        OutputFile.write(f'''\n''')

