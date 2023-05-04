import numpy as np
import pyvisa
import serial
import socket
import json
import time
import sys

rm = pyvisa.ResourceManager()
def numtostr(mystr):`
    return '%20.15e' % mystr

#SR lock-in amplifier
class SR830:
    def __init__(self):
        self._gpib = rm.open_resource('GPIB0::8::INSTR')
    def __del__(self):
        self._gpib.close()
    def mag(self):
        self._gpib.write('OUTP? 3')
        x = self._gpib.read()
        v = np.float_(np.fromstring(x,dtype=float,sep=' '))
        return v