# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument

### Keysight E5071C vector network analyzer
class E5071C:
    def __init__(self, addr):
        self._gpib = rm.open_resource('TCPIP::'+str(addr)+'::INSTR')
    def __del__(self):
        self._gpib.close()
    def idn(self):
        return(self._gpib.query('*IDN?'))
    def save_s2p(self, filename): #save s2p file
        self._gpib.write(':MMEM:STOR \"test.s2p\"')
        self._gpib.timeout = 120000
        data = self._gpib.query(':MMEM:TRAN? \"test.s2p\"') #\"'+filename+'\"')
        self._gpib.timeout = 120000
        data = data[data.index('!'):]
        data = data.replace('\r','')
        with open(filename,'w') as f:
            f.write(data[data.index('!'):])
    def trig(self):
        self._gpib.write(':INIT')
    def avecle(self):
        self._gpib.write('SENS:AVER:CLE')
    def Pow(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SOUR:POW ' + p)
    def cent(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:FREQ:CENT ' + p)
    def marker_max(self):
        self._gpib.write('CALC:MARK:FUNC:TYPE MAX')
        self._gpib.write('CALC:MARK:FUNC:EXEC')
        self._gpib.write('CALC:MARK:SET CENT')
        
    def marker_min(self):
        self._gpib.write('CALC:MARK:FUNC:TYPE MIN')
        self._gpib.write('CALC:MARK:FUNC:EXEC')
        self._gpib.write('CALC:MARK:SET CENT')
        
    def span(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:FREQ:SPAN ' + p)
    def start(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:FREQ:START ' + p)
    def stop(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:FREQ:STOP ' + p)
    def IFBW(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:BWID ' + p )
    def PowerOff(self):
        self._gpib.write('OUTP 0')
        return
    def PowerOn(self):
        self._gpib.write('OUTP 1')
        return
    def SetPoints(self, x):                                 # Not available in PNA-X (Keisght familiy)
        mystr = '%d' % x
        mystr = 'SENS:SWE:POIN ' + mystr
        self._gpib.write(mystr)
    def maxfreq(self):
        data = self._gpib.query('CALC:MARK:X?')
        return(data)