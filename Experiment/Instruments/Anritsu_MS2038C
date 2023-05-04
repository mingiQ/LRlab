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

## Anritsu MS2038C VNA master  ''' TCPIP[board]::host address[::LAN device name][::INSTR] '''
## TCP/IP = 192.168.0.105
## Tested in / RFstation pc / CF room pc /
class MS2038:
    def __init__(self, addr):
        self._gpib = rm.open_resource('TCPIP::'+str(addr)+'::INSTR')
    def mode(self, mod): 
        self._gpib.write(':INST[:SEL] ' + str(mod))  # SPA : spectrum analyzer, MWVNA: vector network analyzer
        #self._gpib.baud_rate = 57600
        return(self._gpib.query(':INST?'))
    def __del__(self):
        self._gpib.close()
    def idn(self):
        return(self._gpib.query('*IDN?'))
    def avgclear(self):
        self._gpib.write(':SENS:AVER:CLE')
    def avgcount(self, p):
        self._gpib.write(':SENS:AVER:COUN ' + str(p))
    def start(self, p):
        self._gpib.write(':SENS:FREQ:STAR ' + str(p))
        self._gpib.baud_rate = 57600
    def stop(self, p):
        self._gpib.write(':SENS:FREQ:STOP ' + str(p))
        self._gpib.baud_rate = 57600
    def startq(self):
        st = self._gpib.query(':SENS:FREQ:STAR?')
        print(st)
    def sweep(self, sweeptype):
        if sweeptype =='single':
            self._gpib.write(':INITiate:HOLD OFF')
            self._gpib.write(':SENSe:SWEep:TYPE SINGle')
        elif sweeptype == 'cont':
            self._gpib.write(':INITiate:HOLD OFF')
            self._gpib.write(':SENSe:SWEep:TYPE CONTinuous')  # defalut val
    def IFBW(self, p):
        self._gpib.write(':SENS:SWE:IFBW ' + str(p))  # must be one of 100000|50000|20000|10000|5000|2000|1000|500|200|100|50|20|10
        self._gpib.baud_rate = 57600
    def points(self, p):
        self._gpib.write(':SENSe:SWE:POINt ' + str(p))
        self._gpib.baud_rate = 57600
    def minimum(self):
        mini = self._gpib.write(':CALCulate:MARKer 1 :MINimum')
        print(mini)
    def save_s2p(self, filename): #save s2p file
        self._gpib.write(':MMEMory:MSIS INTernal')                  # VERY important to define temp data storage first!
        self._gpib.timeout = 120000                                 # Anritsu is really lazy.. give it sufficient time to respond...:( 
        self._gpib.write(':MMEM:STOR:TRAC 4, \"tempdata\"')
        self._gpib.baud_rate = 57600
        time.sleep(1)
        self._gpib.timeout = 120000
        data = self._gpib.query(':MMEM:DATA? \"tempdata.s2p\"')
        self._gpib.baud_rate = 57600
        self._gpib.timeout = 120000
        time.sleep(2)
        data = data.replace('\r','')
        data = data[data.index('!'):]
        time.sleep(1)
        with open(filename,'w') as f:
            f.write(data[data.index('!'):])
            
    def save_trace(self,filename):
        temp = time.strftime("%d%b%Y%H%M%S", time.localtime())
        self._gpib.write(':MMEMory:MSIS INTernal')
        self._gpib.timeout = 120000     
        self._gpib.write(':MMEM:STOR:TRAC 0,\"'+str(temp)+'\"')
        self._gpib.baud_rate = 57600
        self._gpib.timeout = 120000   
        data = self._gpib.query(':MMEMory:DATA? \"'+str(temp)+'.spa\"')
        self._gpib.baud_rate = 57600
        self._gpib.timeout = 120000
        data = data.replace('\r','')
        data = data.replace('MHz', '')
        data = data.replace('=', ',')
        data = data[data.index('\nP_'):]
        data = data[:data.index('#')]
        with open(filename,'w') as f:
            f.write(data.replace('P_', ''))
            
    def transfer(self, filename, path): # transfer internal memory to pc
        data = self._gpib.query(':MMEM:DATA? \"'+ filename +'\"')
        data = data.replace('\r','')
        data = data[data.index('!'):]
        with open(path+filename,'w') as f:
            f.write(data[data.index('!'):])