#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, socket, pickle

# To terminate a DPA-3 Server
# usage: python Terminator.py

ServerHost = socket.gethostname()
ServerPort = 31079 # A prime number brings good luck
MaxBufferSize = 1024
Terminating_Code = -1
ClientSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
ClientSocket.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, MaxBufferSize)
ClientSocket.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, MaxBufferSize)
Sent = ClientSocket.sendto(pickle.dumps(Terminating_Code), (ServerHost, ServerPort))
ClientSocket.close()
