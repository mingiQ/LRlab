# In[0] : necessary packages
import sys
import time
import numpy as np

sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument

### EG&G 5210 Lock-in Amplifier
class EGG5210(VisaInstrument):
    
    def __init__(self, name='e5210', address='GPIB0::15::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=0.05
        self.timeout = 16000
        self.term_char=''
 
    
    def get_idn(self):
        self.instrument.timeout = self.timeout
        return(self.instrument.query("ID"))
    
    def set_voltage(self, v):
        self.instrument.write('OA '+str(v))  # mV
        
    def set_frequency(self, f):
        decimal = int(np.around(np.log10(f),0))
        unit = 4-decimal
        F = int(f*10**(unit))
        BW = str(decimal)
        amp = str(F)
        self.instrument.write('OF '+amp+' '+BW)
        
    def set_sensitivity(self, sen='1V'):
        allowed = ['100nV', '300nV', '1uV', '3uV', '10uV', '30uV', '100uV', '300uV', '1mV', '3mV', '10mV', '30mV', '100mV', '300mV', '1V', '3V ']
        if sen not in allowed:
            print('sensitivity need to be one of ' + ','.join(allowed))
        
        N = allowed.index(sen)
        self.instrument.write('SEN [%d]' % N)
    
    def get_display1(self):
    
        self.instrument.write('D1 [5]')
        ds1 = self.instrument.read()
        return ds1
        
    def meas_X(self):
        #self.instrument.timeout = self.timeout
        self.instrument.write('X')
        x_v = self.instrument.read()
        v = np.float_(np.fromstring(x_v ,dtype=float,sep=' '))    # mV
        return v
    def meas_Y(self):
        self.instrument.write('Y')
        y_v = self.instrument.read()
        v = np.float_(np.fromstring(y_v ,dtype=float,sep=' '))  # mV
        return v  
    def meas_XY(self):
        self.instrument.write('XY')
        xy = self.instrument.read()
        v = np.float_(np.fromstring(xy ,dtype=float,sep=' '))  # mV
        return v
    def meas_mag(self):
        self.instrument.write('MAG')
        mag = self.instrument.read()
        v = np.float_(np.fromstring(mag ,dtype=float,sep=' '))  # mV
        return v 
    def meas_phase(self):
        self.instrument.write('PHA')
        phase = self.instrument.read()
        ph = np.float_(np.fromstring(phase ,dtype=float,sep=' '))  # millidegree
        return ph/1e3
    
    def res_convert(sample, sens, curr):
        full_scale = 10000
        
        meas = sample.meas_X()
        voltage = meas*sens/full_scale  # nV
        resistance = voltage/curr       #Ohm
    
        return voltage, resistance
# =============================================================================
# # In[]
# LIA = EGG5210()
# print(LIA.get_idn())
# 
# =============================================================================
