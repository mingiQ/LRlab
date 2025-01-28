# In[0] : necessary packages
import sys
import time
import numpy as np
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument


### Keithley digital multimeter 2612
class Keithley_2612:
    
    def __init__(self, name='K2612', address='GPIB0::26::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=0.05
        self.timeout = 160000
        self.term_char=''
        
    def get_idn(self):
        return(self.instrument.query("*IDN?"))
    
    def set_voltage(self, volt, ch='a'):
        self.instrument.write(f'smu{ch}.source.levelv = {volt}')
        
    def set_current(self, curr, ch='a'):
        self.instrument.write(f'smu{ch}.source.leveli = {curr}')
    
    def set_source(self, state):
        if state:
            s = "ON"
        else:
            s ="OFF"
        
        self.instrument.write(f'smua.source.output =smua.OUTPUT_{s}')
        
    def get_current(self, delay):
        self.instrument.write('smua.nvbuffer1.clear()')
        self.instrument.write('smua.measure.i(smua.nvbuffer1)')
        time.sleep(delay)
        self.instrument.write('printbuffer(1,1, smua.nvbuffer1.readings)')
        x = self.instrument.read()
        return np.float_(np.fromstring(x,dtype=float,sep=' '))
    
    def get_voltage(self, delay):
        self.instrument.write('smua.nvbuffer1.clear()')
        self.instrument.write('smua.measure.v(smua.nvbuffer1)')
        time.sleep(delay)
        self.instrument.write('printbuffer(1,1, smua.nvbuffer1.readings)')
        x = self.instrument.read()
        return np.float_(np.fromstring(x,dtype=float,sep=' '))
    
    def linearsweep_I(self, starti, stopi, stime, points):
        '''
        SweepILinMeasureV(smua, 1e-3, 10e-3, 0.1, 10)
        --Linear staircase sweep, Channel A, 1mA to 10mA, 0.1 second delay, 10 points.
        '''
        self.instrument.write("smua.reset()")
        self.instrument.write(f"SweepILinMeasureV(smua, {starti}, {stopi}, {stime}, {points})")
        self.instrument.write("waitcomplete()")
        self.instrument.timeout = 160000
        self.instrument.write(f"printbuffer(1, {points}, smua.nvbuffer1.readings)")
        x = self.instrument.read()
        
        return np.fromstring(x.split('\n')[0], dtype=float, sep=',')
    
    def linearsweep_V(self, startv, stopv, stime, points):
        
        '''
        SweepVLinMeasureI(smu, startv, stopv, stime, points):
            Smu: smua for channel A or smub for channel B.
            Start voltage value in volts.
            Stop voltage value in volts.
            Settling time (source-measure delay in seconds).
            Number of points (≥2).

        '''
        self.instrument.write("smua.reset()")
        self.instrument.write(f"SweepVLinMeasureI(smu, {startv}, {stopv}, {stime}, {points})")
        self.instrument.write("waitcomplete()")
        self.instrument.timeout = 160000
        self.instrument.write(f"printbuffer(1, {points}, smua.nvbuffer1.readings)")
        x = self.instrument.read()
        
        return np.fromstring(x.split('\n')[0], dtype=float, sep=',')
    
    def listsweep_I(self, ilist, stime, points):
        '''
        SweepIListMeasureV(smu, ilist, stime, points)
        
        Define current list sweep:
        Smu: smua for channel A or smub for channel B.
        List of current values in amps.
        Settling time (source-measure delay in seconds).
        Number of points (≥2)

        '''
        
        self.instrument.write("smua.reset()")
        
        Ilist = {i for i in ilist}
        points = len(ilist)
        
        self.instrument.write(f"SweepIListMeasureV(smu, {Ilist}, {stime}, {points})")
        self.instrument.write("waitcomplete()")
        self.instrument.timeout = 160000
        self.instrument.write(f"printbuffer(1, {points}, smua.nvbuffer1.readings)")
        x = self.instrument.read()
        
        return np.fromstring(x.split('\n')[0], dtype=float, sep=',')
        
        
    def listsweep_V(self, vlist, stime):
        '''
        SweepVListMeasureI(smu, vlist, stime, points)
        
        Define current list sweep:
        Smu: smua for channel A or smub for channel B.
        List of voltage values in volts.
        Settling time (source-measure delay in seconds).
        Number of points (≥2)

        '''
        
        self.instrument.write("smua.reset()")
        self._DigitaIO_toggle()
        
        Vlist = {v for v in vlist}
        points = len(vlist)
        
        self.instrument.write(f"SweepVListMeasureI(smu, {Vlist}, {stime}, {points})")
        self.instrument.write("waitcomplete()")
        self.instrument.timeout = 160000
        self.instrument.write(f"printbuffer(1, {points}, smua.nvbuffer1.readings)")
        x = self.instrument.read()
        
        return np.fromstring(x.split('\n')[0], dtype=float, sep=',')
    
    def _DigitaIO_toggle(self):
        self.instrument.write("digio.writeibt(1,1)")
        self.instrument.write("digio.writebit(1,0)")
        
# =============================================================================
#     
# # In[]
# #np.fromstring(x,dtype=float,sep=',')
# dmm = Keithley_2612()
# 
# # In[]
# 
# print(dmm.get_idn())
# 
# # In[]
# 
# 
# vlist = np.array([0, 0.1, 0.2, 0.1, 0, -0.1, -0.2, -0.1,0])
# 
# vs  = dmm.listsweep_V(vlist=vlist, stime=0.1)
# 
# 
# # In[]
# 
# V_max = 1
# V_min = -1
# V_step = 0.1
# 
# vlist1 = np.arange(0, V_max, V_step)
# vlist2 = np.arange(V_max, V_min, -V_step)
# vlist3 = np.arange(V_min, V_step, V_step)
# 
# vlist = np.around(np.concatenate((vlist1, vlist2, vlist3)),2)
# 
# vs  = dmm.listsweep_V(vlist=vlist, stime=0.1)
# 
# # In[]
# 
# vs = dmm.linearsweep_V(startv=0, stopv=0.1, stime=0.1, points=10)
# 
# print(vs)
# 
# # In[]
# 
# vs=np.fromstring(vs.split('\n')[0], dtype=float, sep=',')
# 
# # =============================================================================
# # # In[]        
# #     def __init__(self, addr):
# #         self._gpib = rm.open_resource('GPIB0::'+str(addr)+'::INSTR')
# #     def __del__(self):
# #         self._gpib.close()
# #     def setvolt(self,v): # numpy float v
# #         s = np.array2string(v,precision=6)
# #         self._gpib.write('smua.source.levelv = '+s)
# #     def volt_on(self):
# #         self._gpib.write('smua.source.output =smua.OUTPUT_ON')
# #     def volt_off(self):
# #         self._gpib.write('smua.source.output =smua.OUTPUT_OFF')
# #     def meas_i(self, t):
# #         self._gpib.write('smua.nvbuffer1.clear()')
# #         self._gpib.write('smua.measure.i(smua.nvbuffer1)')
# #         time.sleep(t)
# #         self._gpib.write('printbuffer(1,1, smua.nvbuffer1.readings)')
# #         x = self._gpib.read()
# #         v = np.float_(np.fromstring(x,dtype=float,sep=' '))
# #         return v    
# # 
# # =============================================================================
# =============================================================================
