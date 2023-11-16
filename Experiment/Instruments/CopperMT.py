# In[0] : necessary packages
import sys
import time
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.Instruments.InstrumentType import VisaInstrument


## Copper Mountain VNA

class CMT(VisaInstrument):
    def __init__(self, name='CopperMountain', address='TCPIP0::127.0.0.1::5025::SOCKET', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=1
        self.recv_length=65536
        self.term_char=''
        
        
# In[]
        
        RM =  pyvisa.ResourceManager('C:/Windows/System32/visa32.dll')
        self._gpib = RM.open_resource('TCPIP0::127.0.0.1::5025::SOCKET')
    
    def start(self, freq):
        p = numtostr(freq)
        self._gpib.write(':SENS1:FREQ:STAR '+p+' \n')
    def stop(self, freq):
        p = numtostr(freq)
        self._gpib.write(':SENS1:FREQ:STOP '+p+' \n')
    def points(self, pnts):
        p = numtostr(pnts)
        self._gpib.write(':SENS1:SWE:POIN '+p+' \n')
    def avgcount(self, avg):
        p = numtostr(avg)
        self._gpib.write('SENS1:AVER:COUN' +p+' \n')
    def avgclear(self):
        self._gpib.write('SENS1:AVER:CLE \n')
    def IFBW(self, freq):
        p = numtostr(freq)
        self._gpib.write(':SENS1:BWID '+p+' \n')
    def power(self, po):
        p = numtostr(po)
        self._gpib.write(':SOUR:POW '+p+' dBm \n')
    def sweep(self, mode):
        if mode == 'single':
            self._gpib.write(':INIT  \n')
            self._gpib.write(':TRIG:SOUR BUS \n')
            self._gpib.write(':TRIG:SEQ:SING \n')
    def save_s2p(self, filename):
        self._gpib.write(':MMEM:STOR:SNP \"test2.s2p\" \n')
        self._gpib.timeout = 120000
        with open("C:/Program Files/S2VNA/FixtureSim/test2.s2p", 'r') as source, open(filename, 'w') as target:
            for line in source:
                target.write(line)
        