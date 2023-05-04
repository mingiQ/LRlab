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

### Keysight 34461A digital voltmeter
class DMM:
    def __init__(self, addr):
        self._gpib = rm.open_resource('GPIB0::'+str(addr)+'::INSTR')
    def __del__(self):
        self._gpib.close()
    def volt(self):
        self._gpib.write('READ?')
        x = self._gpib.read()
        v = np.float_(np.fromstring(x,dtype=float,sep=' '))
        return v
    def null(self, s):
        self._gpib.write('SENS:VOLT:NULL ' + str(s))