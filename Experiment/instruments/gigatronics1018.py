# In[0] : necessary packages
import sys
import time
import numpy as np

sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument

### EG&G 5210 Lock-in Amplifier
class GT1018(VisaInstrument):
    def __init__(self, name='gt5210', address='GPIB0::8::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=0.05
        self.timeout = 160000
        self.term_char=''
        #self.prefix = '$'
        
    def flush(self):
        self.instrument.timeout = 160000
        self.instrument.write('SEND NUL')
    
    def get_idn(self):
        self.instrument.write('SEND STATUS')
        time.sleep(self.query_sleep)
        msg = self.instrument.read()
        return print(f"{msg} Gigatronics 1018")
        
    def set_frequency(self, freq):
        '''
        Parameters
        ----------
        freq : MHz

        set output frequency
        -------

        '''
        self.instrument.timeout = 160000
        self.instrument.write('GEN FIXED')
        self.instrument.write(f'FA {freq}')
        
    def get_frequency(self):
        
        
        self.instrument.write('SEND FREQ')
        time.sleep(self.query_sleep)
        msg = self.instrument.read()
        return float(msg.split('F OUT')[1].split('\r\n')[0])
        
    def set_power(self, power):
        '''
        Parameters
        ----------
        power : dBm

        set output power
        -------

        '''
        self.instrument.timeout = 160000
        self.instrument.write(f'LEVEL {power}')
        
    def get_power(self):
        
        self.instrument.write('SEND POWER')
        time.sleep(self.query_sleep)
        msg = msg = self.instrument.read()
        return float(msg.split('P INT')[1].split('\r\n')[0])
    
    def set_sweep(self, f_start, f_stop, power, f_step=1, rate=10, swpmode='single'):
        '''
        sweep time: 10ms - 100s per step, 10 options
        step: 1MHz, 10MHz, 100MHz
        
        rate: 1 : 1.0 MHz/s
        1.0 steps/s
        rate: 2 : 2.4 MHz/s
        2.4 steps/s
        rate: 3 : 5.0 MHz/s
        5.0 steps/s
        rate: 4 : 10.0 MHz/s
        10.0 steps/s
        rate: 5 : 20.2 MHz/s
        20.2 steps/s
        rate: 6 : 57.8 MHz/s
        57.8 steps/s
        rate: 7 : 89.0 MHz/s
        89.0 steps/s
        rate: 8 : 178.2 MHz/s
        178.2 steps/s
        rate: 9 : 1200.0 MHz/s
        1200.0 steps/s
        rate: 10 : 1200.0 MHz/s
        1200.0 steps/s

        Parameters
        ----------
        f_start : [MHz]
        f_stop : [MHz]
        f_step : [MHz] [front pannel: sweep increment]
        rate : 1 (slowest) - 10 (fastest) [front pannel: sweep time]
        Returns
        -------
        None.

        '''
        swprate = {1:'A', 2:'B', 3:'C', 4:'D', 5:'E', 6:'F', 7:'G', 8:'H', 9:'I', 10:'J'}
        mode = {'single':'ONCE', 'step':'STEP', 'auto':'AUTO'}
        
        
        
        self.instrument.write('GEN USWP')
        self.instrument.write(f'FA {f_start}')
        self.instrument.write(f'FB {f_stop}')
        self.instrument.write(f'FC {f_step}')
        self.instrument.write(f'LEVEL {power}')
        self.instrument.write(f'SWPRATE {swprate[rate]}')
        self.instrument.write(f'SWEEP {mode[swpmode]}')
        
        
        
        
    def set_sweep_terminate(self):
        self.instrument.write('SWEEP RESET')
        

# =============================================================================
# # In[]
# 
# sg = GT1018()
# 
# # In[]    
# 
# sg.set_frequency(3000)
# sg.set_power(-90)
# 
# time.sleep(0.1)
# p = sg.get_power()
# f = sg.get_frequency()
# 
# print(f"freq={f}MHz, {p}dBm")
# 
# # In[]
# 
# import pyvisa
# 
# RM =  pyvisa.ResourceManager('C:/Windows/System32/visa32.dll')
# gt= RM.open_resource('GPIB0::8::INSTR')       
# 
# 
# # In[]
# import time
# 
# t = 0.05
# 
# gt.timeout = 160000
# gt.write('FA 1000')
# #gt.write('SEND NUL')
# #gt.flush()
# gt.timeout = 160000
# gt.write('SEND FREQ')
# #gt.write('SEND FREQ')
# time.sleep(t)
# #gt.read_raw()
# a = gt.read()
# print(a)
# 
# # In[]
# 
# gt.write('SEND STATUS')
# #gt.write('SEND FREQ')
# time.sleep(t)
# #gt.read_raw()
# a = gt.read()
# print(a)
# 
# # In[]
# gt.timeout = 160000
# gt.write('LEVEL -100')
# #time.sleep(1)
# #gt.clear()
# gt.write('SEND POWER')
# #gt.write('SEND POWER')
# time.sleep(t)
# #gt.read_raw()
# a = gt.read()
# print(a)
# 
# # In[]
# gt.clear()
# 
# 
# # In[rate cal]
# 
# #rate = 10
# step = 1
# t_wait = 5
#     
# for rate in [8,9,10]: #np.linspace(1,10,10):
#     
#     rate = int(rate)
#     sg.set_sweep(f_start=1000, f_stop=17000, power=-20, f_step=step, rate=rate, swpmode='single')
#     
#     time.sleep(1)
#     
#     t0 = time.time()
#     f1 = sg.get_frequency()
#     t1 = time.time()
#     
#     time.sleep(t_wait)
#     
#     f2 = sg.get_frequency()
#     t2 = time.time()
#     
#     sg.set_sweep_terminate()
#     
#     ft = (f2-f1)/(t2-t1)
#     print(f'rate: {rate} : {ft} MHz/s')
#     print(f'{(t1-t0)*1e3}ms : gpib response time')
#     print(f'{ft/step} steps/s')
# 
# # In[]
# from binascii import hexlify
# 
# s = 'M'.encode()
# prefix = '$'    
# print(prefix+hexlify(s, ',' ).decode())
# 
#   
# # In[]
# 
# s = 'MESSAGE'
# print('%d'*len(s) % tuple(map(hex,map(ord, s))))
# 
# 
# 
# 
# 
# 
# =============================================================================
