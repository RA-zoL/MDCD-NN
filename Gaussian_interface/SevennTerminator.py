#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import socket, pickle

TerminatingCode = -1
Host = socket.gethostname()
Port = 31079 # A prime number brings good luck
MaxBufferSize = 1024
ClientSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
ClientSocket.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, MaxBufferSize)
ClientSocket.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, MaxBufferSize)
Sent = ClientSocket.sendto(pickle.dumps(TerminatingCode), (Host, Port))
Data, ServerAddress = ClientSocket.recvfrom(MaxBufferSize)
ClientSocket.close()
