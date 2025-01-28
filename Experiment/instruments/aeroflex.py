# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import  VisaInstrument

import time
import numpy as np
import glob
import os.path


## Weinschel programmable step attenuator 


class Aeroflex(VisaInstrument):


    def __init__(self, name="aeroflex", address='GPIB::10::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled=enabled)
        self.query_sleep = 0.05
        self.timeout = 10000
        self.term_char=''
        self.instrument.read_termination='\n'
        

    def get_idn(self):
        return self.instrument.query('*IDN?')
    

    def get_query_sleep(self):
        return self.query_sleep
    
    def set_ch(self, ch, val):
        self.instrument.write(f'ATTN CH{ch} {val}')
        
    def read_ch(self, ch):
        return (self.instrument.query(f'ATTN? CH{ch}'))
    


    def get_settings(self):
        settings = {"ch1": self.read_ch(1), 
                    "ch2": self.read_ch(2), 
                    }
        return settings
    
    

# =============================================================================
# # In[]
# 
# if __name__ == '__main__':
#     na = PNA("N5242A", address="192.168.0.102")
#     print(na.get_id())
#     
#     
#     
# # In[]
# import pyvisa
# # 
# RM =  pyvisa.ResourceManager()
# na = RM.open_resource('TCPIP0::192.168.0.102::INSTR')   
# 
# # In[]
# 
# na = RM.open_resource('TCPIP0::K-N5242B-22393::hislip0,4880::INSTR')   
# 
# 
# # In[]
# #na.timeout=160000
# idn = na.query('*IDN?')
# print(idn)
# 
# na.write('CALC:PAR:SEL CH1_S11_1')
# 
# # In[]
# 
# na.write("CALC1:FORM MLOG")
# # In[]
# na.timeout=16000
# data = na.query("CALC1:FORM?")
# 
# # In[]
# na.timeout=100
# fffd = na.query("CALC1:DATA? FDATA")
# 
# # In[]
# 
# 
# data = na.query(":SENS1:FREQ:STOP?")
# 
# 
# =============================================================================
