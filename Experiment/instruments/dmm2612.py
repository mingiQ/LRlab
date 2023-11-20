# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument


### Keithley digital multimeter 2612
class Keithley_2612:
    def __init__(self, addr):
        self._gpib = rm.open_resource('GPIB0::'+str(addr)+'::INSTR')
    def __del__(self):
        self._gpib.close()
    def setvolt(self,v): # numpy float v
        s = np.array2string(v,precision=6)
        self._gpib.write('smua.source.levelv = '+s)
    def volt_on(self):
        self._gpib.write('smua.source.output =smua.OUTPUT_ON')
    def volt_off(self):
        self._gpib.write('smua.source.output =smua.OUTPUT_OFF')
    def meas_i(self, t):
        self._gpib.write('smua.nvbuffer1.clear()')
        self._gpib.write('smua.measure.i(smua.nvbuffer1)')
        time.sleep(t)
        self._gpib.write('printbuffer(1,1, smua.nvbuffer1.readings)')
        x = self._gpib.read()
        v = np.float_(np.fromstring(x,dtype=float,sep=' '))
        return v    
