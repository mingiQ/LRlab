# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import SocketInstrument, VisaInstrument

import time
import numpy as np
import glob
import os.path

#
# ======================================================
#
# reference the programming manual here: http://na.support.keysight.com/pna/help/index.html?id=1000001808-1:epsg:man
# and specifically here: New Programming Commands:
# http://na.support.keysight.com/pna/help/latest/help.htm
# help: 1 800 829-4444
#
# ======================================================
#

def polar2mag(xs, ys):
    return 20*np.log10(np.sqrt(xs ** 2 + ys ** 2)), np.arctan2(ys, xs) * 180/np.pi

class E4421B(VisaInstrument):
    """
    The interface to the Agilent E4421B RF Generator, implemented on top of 
    :py:class:`~slab.instruments.instrumenttypes.SocketInstrument`
    """
    def __init__(self, name='E4421B', address='GPIB0::24::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=0.05
        self.timeout = 160000
    
    def get_id(self):
        """Get Instrument ID String"""
        return self.instrument.query('*IDN?').strip()
    
    def set_output(self,state=True):
        """Set Output State On/Off"""
        if state: self.instrument.write(':OUTPUT:STATE ON')
        else:     self.instrument.write(':OUTPUT:STATE OFF')
        
    def get_output(self):
        """Query Output State"""
        return int(self.instrument.query(':OUTPUT?')) == 1
    
    def set_mod(self,state=True):
        """Set Modulation State On/Off"""
        if state: self.instrument.write(':OUTPUT:MOD:STATE ON')
        else:     self.instrument.write(':OUTPUT:MOD:STATE OFF')
        
    def get_mod(self):
        """Query Modulation State"""
        return bool(self.instrument.query(':OUTPUT:MOD?'))
            
    def set_frequency(self, frequency):
        """Set CW Frequency in Hz"""
        self.instrument.write(':FREQUENCY %f' % frequency)
    
    def get_frequency(self):
        """Query CW Frequency"""
        return float(self.instrument.query(':FREQUENCY?'))

    def set_sweep(self,start,stop,numpts,dwell):
        """Sets up frequency sweeping parameters"""
        self.instrument.write(':FREQUENCY:START %f;:FREQUENCY:STOP %f; :SWEEP:POINTS %f; :SWEEP:DWELL %f' % (start,stop,numpts,dwell))

    def get_sweep(self):
        """Gets current frequency sweeping parameters"""
        return [float(s) for s in (self.instrument.query(':FREQUENCY:START?'),self.instrument.query(':FREQUENCY:STOP?'), self.instrument.query(':SWEEP:POINTS?'),self.instrument.query(':SWEEP:DWELL?'))]

    def set_sweep_mode(self,enabled=True):
        """Set's up source for sweep mode"""
        if enabled:
            self.instrument.write(':LIST:TYPE STEP; :LIST:TRIG:SOURCE EXT; :FREQuency:MODE SWEEP')
        else:
            self.instrument.write(':FREQ:MODE CW; :LIST:TRIG:SOURCE CONT')

    def set_cw_mode(self):
        """Set generator into CW mode"""
        self.instrument.write(':FREQ:MODE CW; :LIST:TRIG:SOURCE CONT')

    def set_phase(self,phase):
        """Set signal Phase in radians"""
        self.instrument.write(':PHASE %f' % phase)
    
    def get_phase(self):
        """Query signal phase in radians"""
        return float(self.instrument.query(':PHASE?'))
        
    def set_power(self,power):
        """Set CW power in dBm"""
        self.instrument.write(':POWER %f' % power)
        
    def get_power(self):
        return float(self.instrument.query(':POWER?'))

    def get_settings(self):
        settings=SocketInstrument.get_settings(self)
        settings['frequency']=self.get_frequency()
        settings['power']=self.get_power()
        settings['phase']=self.get_phase()
        settings['mod']=self.get_mod()
        settings['output']=self.get_output()
        settings['id']=self.get_id()
        return settings
        
    def get_settled(self):
        """Get source settled state"""
        return bool(self.instrument.query(':OUTPut:SETTled?'))
        
    def wait_to_settle(self, timeout=1):
        """Block until source settled"""
        start=time.time()
        while not self.get_settled() and time.time()-start<timeout: pass
    
    def set_internal_pulse(self,pulse_time=10e-6):
        self.instrument.write(':SOUR:PULM:SOUR:INT FRUN')
        self.instrument.write(':SOUR:PULM:INT:PERIOD %f S' %(pulse_time))
        self.instrument.write(':SOUR:PULM:INT:PWIDTH %f S' %(pulse_time))
        self.instrument.write(":SOUR:PULM:STAT ON")
        self.set_mod()
    
    def set_ext_pulse(self,mod=True):
        self.instrument.write(':SOUR:PULM:SOUR EXT')
        if mod:
            self.instrument.write(":SOUR:PULM:STAT ON")
            self.set_mod(True)
        else:
            self.instrument.write(":SOUR:PULM:STAT OFF")
            self.set_mod(False)


# =============================================================================
# # In[]
# 
# 
# sg = E4421B()
# print(sg.get_id())
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
# 
# =============================================================================
