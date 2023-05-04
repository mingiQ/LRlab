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

## DAC LR
class DAC_CF:
    def __init__(self, addr):       # COM4
        self.DAC = serial.Serial(port=str(addr), baudrate=115200, timeout=.1)
    def idn(self):
        self.DAC.write(bytes('*IDN?', 'utf-8'))
        data = self.DAC.readline().decode('utf-8').rstrip()
        print(data)
    def voltsweep(self, dac_CF, volt):
        self.DAC.write(bytes("WF "+ str(dac_CF) +" "+ str(volt) +" "+ str(volt), 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().decode('utf-8').rstrip()
        print(data)