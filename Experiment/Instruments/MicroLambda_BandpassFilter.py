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

## MicroLambda bandpass filter --> need to fix!
class MLBPF:
    def __init__(self, addr):
        self.arduino = serial.Serial(port=str(addr), baudrate=115200, timeout=0.1)
    def write_read(self, x):
        self.arduino.write(bytes(x, 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().decode('utf-8').rstrip()
        print(data)
    def idn(self):
        self.arduino.write(bytes('?', 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().decode('utf-8').rstrip()
        print(data)
    def freq(self, f):
        self.arduino.write(bytes('FH' + ' ' + str(f), 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().decode('utf-8').rstrip()
        print(data)
    def quit_bpf(self):
        self.arduino.close()