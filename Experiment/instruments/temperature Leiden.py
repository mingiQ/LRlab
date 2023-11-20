import numpy as np
import pyvisa
import serial
import socket
import json
import time
import sys


### Leiden temperature monitor
class Thermo_leiden:
    def __init__(self):
        self.HOST = '10.164.26.139'
        self.PORT = 8888
        self.soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.soc.connect((self.HOST, self.PORT))
        self.repl = self.soc.recv(1024)
        print(self.repl.decode())
        
    def readT(self, addr):
        self.cmd = str(addr)                #3K	Still	CP	MC TT-1701	3k-pr	Still-Pr/Crnx/PT	CP-50mK-Pr	MC-pr	CMN 089 MC	CMN 060 pr
        self.soc.send(self.cmd.encode())
        repl = self.soc.recv(1024)
        return float(repl)
    
    def Tcal(self, addr):
        self.cmd = 'CAL' + str(addr)
        self.soc.send(self.cmd.encode())
        repl = self.soc.recv(1024)
        return float(repl)
    
    def quit_Tc(self):
        self.cmd = 'quit'
        self.soc.send(self.cmd.encode())
        
    def close_Tc(self):
        self.cmd = 'close'
        self.soc.send(self.cmd.encode())