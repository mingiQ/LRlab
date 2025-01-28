# In[0] : necessary packages
import sys
import time
import numpy as np
import pyvisa
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument


### Keithley sourcemeter 237

class Keithley_237:
    
    def __init__(self, name='K237', address='GPIB0::26::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=0.05
        self.timeout = 160000
        self.term_char=''
        
    def get_idn(self):
        return(self.instrument.query("*IDN?"))
    
    def outputOn(self):
        self.instrument.write('N1X')
        
    def outputOff(self):
        self.instrument.write('N0X')
    
    def set_voltage(self, volt, delay=0):
        if volt == 0:
            v_range = 1
        elif int(np.log10(abs(volt)/0.11)) <= 0:
            v_range = 1
        else:
            v_range = int(np.log10(abs(volt)/0.11))+1

        self.instrument.write(f'F0,0X')
        if v_range > 4:
            v_range = 4
            
        self.instrument.write(f'B{volt},{v_range},{delay}X')
        self.instrument.write('H0X')
        
    def set_current(self, curr, delay=0):
        if curr != 0:
            i_range = int(np.log10(abs(curr)))+10
        elif curr == 0:
            i_range = 1
            
        self.instrument.write(f'F1,0X')
        if i_range > 10:
            i_range = 10
            
        self.instrument.write(f'B{curr},{i_range},{delay}X')
        self.instrument.write('H0X')
    
        
    def get_source(self, source, delay=0):
            
        self.instrument.write(f'G{1},{0},{0}X')

        time.sleep(delay)

        x = self.instrument.read()
        return np.float64(x.split(',B')[0].split(f'NSDC{source}')[1])
    
    def get_meas(self, source, delay=0):
        measure = {'I': 'V', 'V': 'I'}
        self.instrument.write(f'G{4},{0},{0}X')

        time.sleep(delay)

        x = self.instrument.read()
        return np.float64(x.split(',B')[0].split(f'NMDC{measure[source]}')[1])
    
    def flush(self):
        self.instrument.flush(pyvisa.constants.VI_WRITE_BUF_DISCARD)
    
# =============================================================================
# # In[]
# import pyvisa
# rm = pyvisa.ResourceManager()
# #dmm = Keithley_237(address='GPIB0::30::INSTR')
# 
# dmm = rm.open_resource('GPIB0::30::INSTR')
# =============================================================================
