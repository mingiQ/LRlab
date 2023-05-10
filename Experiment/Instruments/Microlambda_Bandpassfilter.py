# In[0] : necessary packages
import numpy as np
import pyvisa
import serial
import socket
import json
import time
import sys

## MicroLambda bandpass filter
class MLBPF:
    '''
        /*    ------ communication commands ---------
     *     serial @ 115200, line end character: \n
     *    "?" or "*IDN?"             ->  devstamp         - equipment ID  
     *    "FF val"                   -> "FF"              - set single dac value (float)
     *    "SF fmin fmax"             -> "SF"              - save fmin and fmax values to EEPROM
     *    "RF"                       -> "RF"              - read fmin and fmax values from EEPROM
     *    
     *    "TS dac dval"              -> "TS"              - comm test: write xAAAAA to dac dval times
     *    "TX pin dval"              -> "TX"              - pin test: write HIGH/LOW to pin dval times
     *    "TT pin"                   -> "TT"              - toggle pin
    */
    '''
    def __init__(self, addr):
        self.arduino = serial.Serial(port=str(addr), baudrate=115200, timeout=None)
        self.arduino.flushInput()
        self.arduino.flushOutput()
    def write_read(self, x):
        self.arduino.write(bytes(x, 'utf-8'))
        time.sleep(0.05)
        #data = self.arduino.readline().decode('utf-8').rstrip()
        data = self.arduino.readline().rstrip()
        print(data)
    def idn(self):
        self.arduino.write(bytes('?', 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().rstrip()
        print(data)
    def freq(self, f):
        self.arduino.write(bytes('FF' + ' ' + str(f), 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().rstrip()
        print(data)
    def debug(self, f):
        self.arduino.write(bytes('FH' + ' ' + str(f), 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().rstrip()
        print(data)
    def quit_bpf(self):
        self.arduino.close()
