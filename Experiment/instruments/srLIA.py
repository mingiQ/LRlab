# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument

#SR lock-in amplifier
class SR830:
    def __init__(self):
        self._gpib = rm.open_resource('GPIB0::8::INSTR')
    def __del__(self):
        self._gpib.close()
    def mag(self):
        self._gpib.write('OUTP? 3')
        x = self._gpib.read()
        v = np.float_(np.fromstring(x,dtype=float,sep=' '))
        return v