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

## AMI430 magnet

class AMI430:
    def __init__(self, addr):       # COM7
        self.magnet = serial.Serial(port=str(addr), baudrate=115200, timeout=.1)