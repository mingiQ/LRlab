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

### EG&G 5210 Lock-in Amplifier
class EGG5210:
    def __init__(self, addr):
        self._gpib = rm.open_resource('GPIB0::'+str(addr)+'::INSTR')
    def __del__(self):
        self._gpib.close()
    def idn(self):
        return(self._gpib.query("ID"))
    def setvolt(self, v):
        self._gpib.write('OA '+str(v))  # mV
    def setfreq(self, f):
        decimal = int(np.around(np.log10(f),0))
        unit = 4-decimal
        F = int(f*10**(unit))
        BW = str(decimal)
        amp = str(F)
        self._gpib.write('OF '+amp+' '+BW)
    def meas_X(self):
        self._gpib.write('X')
        x_v = self._gpib.read()
        v = np.float_(np.fromstring(x_v ,dtype=float,sep=' '))    # mV
        return v
    def meas_Y(self):
        self._gpib.write('Y')
        y_v = self._gpib.read()
        v = np.float_(np.fromstring(y_v ,dtype=float,sep=' '))  # mV
        return v  
    def meas_XY(self):
        self._gpib.write('XY')
        xy = self._gpib.read()
        v = np.float_(np.fromstring(xy ,dtype=float,sep=' '))  # mV
        return v
    def meas_mag(self):
        self._gpib.write('MAG')
        mag = self._gpib.read()
        v = np.float_(np.fromstring(mag ,dtype=float,sep=' '))  # mV
        return v 
    def meas_phase(self):
        self._gpib.write('PHA')
        phase = self._gpib.read()
        ph = np.float_(np.fromstring(phase ,dtype=float,sep=' '))  # millidegree
        return ph/1e3
    
        