import numpy as np
import pyvisa
import serial
import socket
import json
import time
import sys

class AMI430:
    '''
    ​X	​2042 Ohm	​2020 Ohm	​2.0 A/min	​​0.036 T/min	​0.01825_T/A	​	​	​
    ​Y	​2079 Ohm	​2056 Ohm	​2.0 A/min	0.036 T/min	0.01803_T/A	​	​	​
    ​Z	​1572 Ohm	1555 Ohm​	​1.4 A/min	0.165 T/min​	​0.1164_T/A​	​0.05mApp=0.006mTpp	​	​-0.8mA=-0.09mT
    '''
    def __init__(self, addr):       # COM8: Bz COM6: By COM7: Bx
        self.magnet = serial.Serial(port=str(addr), baudrate=115200, timeout=0.1)
        self.magnet.flushInput()
        self.magnet.flushOutput()
    def write_read(self, x):
        self.magnet.write(str.encode(x+'\n'))
        time.sleep(0.05)
        data = self.magnet.readline().rstrip()
        print(data)
    def self_test(self):
        self.magnet.write(str.encode('*TST?\n'))
        time.sleep(0.05)
        data = self.magnet.readline().rstrip()
        print(data)
    def idn(self):
        self.magnet.write(str.encode('*IDN?\n'))
        time.sleep(0.05)
        data = self.magnet.readline().rstrip()
        print(data.decode())
    def debug(self):
        self.magnet.write(str.encode('SYSTem:ERRor?\n'))
        time.sleep(0.05)
        data = self.magnet.readline().rstrip()
        print(data)
    def current_limit(self):
        self.magnet.write(str.encode('CURRent:LIMit?\n'))
        time.sleep(0.05)
        data = self.magnet.readline().rstrip()
        print(data.decode()+'A')
    def unit_set(self, cmd):
        if cmd == 'T':
            com = 1
        elif cmd == 'kG':
            com = 0
        else:
            print('type proper unit')
            
        self.magnet.write(str.encode('CONFigure:FIELD:UNITS ' + str(com)+'\n'))
        self.magnet.write(str.encode('FIELD:UNITS?\n'))
        time.sleep(0.05)
        ind = self.magnet.readline().rstrip()
        unit = ['kG', 'T']
        
        print('unit set to '+ unit[int(ind.decode())-1]) 
        
    def coilconst(self):
        self.magnet.write(str.encode('COILconst? \n'))  #  0.1164 T/A
        time.sleep(0.05)
        cc = self.magnet.readline().rstrip()
        print(cc.decode() + 'T/A')
        
    def set_field(self, field):
        '''
        Require coil constant
        ----------
        field : in Tesla or kG

        Returns
        set b field in tesla or kG and print the set field
        -------
        '''
        self.magnet.write(str.encode('CONFigure:FIELD:TARGet '+str(field)+'\n'))
        self.magnet.write(str.encode('FIELD:TARGet?\n'))
        time.sleep(0.05)
        target = self.magnet.readline().rstrip()
        print("set field: " + target.decode() + 'T')
        
    def set_ramping(self, seg, rate, upper_bound):
        self.magnet.write(str.encode(f'CONFigure:FIELD:TARGet {seg} {rate} {upper_bound}\n'))
        self.magnet.write(str.encode(f'RAMP:RATE:FIELD:{seg}?\n'))
        time.sleep(0.05)
        target = self.magnet.readline().rstrip()
        print(target.decode())
    
    def start_ramping(self):
        ans = input("Ready to ramp?: ")
        if ans == 'y':
            self.magnet.write(str.encode(f'RAMP\n'))
        else:
            pass
        
    def pause(self):
        self.magnet.write(str.encode(f'PAUSE\n'))
        
    def current_field(self):
        self.magnet.write(str.encode(f'FIELD:MAGnet? \n'))
        time.sleep(0.05)
        field = self.magnet.readline().rstrip()
        print(field.decode()+' T')
    
    def current_status(self):
        self.magnet.write(str.encode(f'STATE? \n'))
        time.sleep(0.05)
        status = self.magnet.readline().rstrip()
        ind = int(status.decode())-1
        
        window = ['RAMPING to target field/current',
                  'HOLDING at the target field/current',
                  'PAUSED',
                  'Ramping in MANUAL UP mode',
                  'Ramping in MANUAL DOWN mode',
                  'ZEROING CURRENT (in progress)',
                  'Quench detected',
                  'At ZERO current',
                  'Heating persistent switch',
                  'Cooling persistent switch'
                    ]
        
        print(window[ind])
    
    def quit_magnet(self):
        self.magnet.close()
     