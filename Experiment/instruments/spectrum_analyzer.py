# -*- coding: utf-8 -*-
# In[0]: necesaray packages
"""
Created on Mon Aug 01 21:26:31 2011

@author: Dave

@modified by mingi (11/22/2024)
"""
import sys
import time
import numpy as np
import glob
import os.path
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument


# In[1]: instrument

class E4440(VisaInstrument):
    MAXSWEEPPTS = 1601
    default_port = 5025

    def __init__(self, name="E4440", address='TCPIP::192.168.0.60::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled=enabled)
        self.query_sleep = 0.05
        self.timeout = 10000
        self.term_char=''
        self.instrument.read_termination='\n'

    def get_id(self):
        return self.instrument.query('*IDN?')

    def get_query_sleep(self):
        return self.query_sleep

    def calibrate(self):
        """
        Performs a full alignment
        :return:
        """
        print("Performing full alignment...")
        self.instrument.write(":CAL:ALL")

    #### Frequency setup
    def set_start_frequency(self, freq):
        self.instrument.write(":SENS:FREQ:START %f" % (freq))

    def get_start_frequency(self):
        return float(self.instrument.query(":SENS:FREQ:START?"))

    def set_stop_frequency(self, freq):
        self.instrument.write(":SENS:FREQ:STOP %f" % (freq))

    def get_stop_frequency(self):
        return float(self.instrument.query(":SENS:FREQ:STOP?"))

    def set_center_frequency(self, freq):
        self.instrument.write(":SENS:FREQ:CENTer %f" % (freq))

    def get_center_frequency(self):
        return float(self.instrument.query(":SENS:FREQ:CENTer?"))

    def set_span(self, span):
        return self.instrument.write(":SENS:FREQ:SPAN %f" % (span))

    def get_span(self):
        return float(self.instrument.query(":SENS:FREQ:SPAN?"))

    def set_sweep_points(self, numpts=8192):
        self.instrument.write(":SENSe:SWEep:POINts %f" % (numpts))

    def get_sweep_points(self):
        return float(self.instrument.query(":SENSe:SWEep:POINts?"))

    #### Averaging
    def set_averages(self, averages):
        self.instrument.write(":SENS:AVERage:COUNt %d" % (averages))

    def get_averages(self):
        return int(self.instrument.query(":SENS:AVERage:COUNt?"))

    def set_average_state(self, state=True):
        if state:
            s = "ON"
        else:
            s = "OFF"
        self.instrument.write(":SENS:AVERage:STATe %s" % (s))

    def get_average_state(self):
        return bool(self.instrument.query(":SENS:AVERage:STATe?")[0])

    def clear_averages(self):
        self.instrument.write(":SENS:average:clear")

    def set_resbw(self, bw, auto=None):
        if auto is not None:
            if auto:
                self.instrument.write(":SENS:BANDwidth:RESolution:AUTO ON")
            else:
                self.instrument.write(":SENS:BANDwidth:RESolution:AUTO OFF")
        else:
            self.instrument.write(":SENS:BANDwidth:RESolution %f" % (bw))

    def get_resbw(self):

        return float(self.instrument.query(":SENS:BANDwidth:RESolution?"))

    def set_vidbw(self, bw, auto=None):
        if auto is not None:
            if auto:
                self.instrument.write(":SENS:BANDwidthVIDEO:AUTO ON")
            else:
                self.instrument.write(":SENS:BANDwidth:VIDEO:AUTO OFF")
        else:
            self.instrument.write(":SENS:BANDwidth:VIDEO %f" % (bw))

    def get_vidbw(self):

        return float(self.instrument.query(":SENS:BANDwidth:VIDEO?"))

    def wait_for_completion(self):
        """
        from the e4440 documentation:
        *OPC? query stops new commands from being processed until the current processing is complete.
        """
        return self.instrument.query("*OPC?")

    def restart_measurement(self):
        self.instrument.write(':INITiate:RESTart')

    def set_continuous_sweep_state(self, state=True):
        if state:
            stateString = 'ON'
        else:
            stateString = 'OFF'
        self.instrument.write(':INITiate:CONTinuous ' + stateString)

    def get_continuous_sweep_state(self, state=True):
        if state:
            stateString = 'ON'
        else:
            stateString = 'OFF'
        return bool(self.instrument.query(':INITiate:CONTinuous?' + stateString))

    def trigger_single(self):
        self.instrument.write(':INIT:IMM')

    def set_trigger_source(self, source="IMMEDIATE"):  #IMMEDIATE, LINE, EXTERNAL
        self.instrument.write('TRIG:SOURCE ' + source)

    def get_trigger_source(self):  #INTERNAL, MANUAL, EXTERNAL,BUS
        return self.instrument.query(':TRIGGER:SOURCE?')

        #### File Operations
        # def save_file(self,fname):
        #     self.instrument.write('MMEMORY:STORE:FDATA \"' + fname + '\"')
        #
        # def set_format(self,trace_format='MLOG',channel=1):
        #     """set_format: valid options are
        #     {MLOGarithmic|PHASe|GDELay| SLINear|SLOGarithmic|SCOMplex|SMITh|SADMittance|PLINear|PLOGarithmic|POLar|MLINear|SWR|REAL| IMAGinary|UPHase|PPHase}
        #     """
        #     self.instrument.write(":CALC:FORMAT "+trace_format)
        # def get_format(self,channel=1):
        #     """set_format: valid options are
        #     {MLOGarithmic|PHASe|GDELay| SLINear|SLOGarithmic|SCOMplex|SMITh|SADMittance|PLINear|PLOGarithmic|POLar|MLINear|SWR|REAL| IMAGinary|UPHase|PPHase}
        #     """
        #     return self.instrument.query(":CALC:FORMAT?")

    def read_data(self, channel=1):
        """Read current NWA Data, return fpts,mags,phases"""
        #        self.instrument.write(":CALC1:PAR1:SEL")
        #        self.instrument.write(":INIT1:CONT OFF")
        #        self.instrument.write(":ABOR")
        self.instrument.write(":FORM:DATA ASC")
        #self.instrument.write(":CALC1:DATA:FDAT?")
        #self.instrument.query(":TRAC:DATA? TRACE1")
        data_str =  self.instrument.query(":TRAC:DATA? TRACE1")

# =============================================================================
#         done = False
#         ii = 0
#         while not done:
#             time.sleep(self.query_sleep)
#             ii += 1
#             try:
#                 s = self.instrument.read()
#             except:
#                 print("read %d failed!" % ii)
#             data_str += s
#             done = data_str[-1] == '\n'
# =============================================================================
        #print data_str
        data = np.fromstring(data_str, dtype=float, sep=',')
        # data=data.reshape((-1,2))
        data = data.transpose()
        #self.data=data
        fpts = np.linspace(self.get_start_frequency(), self.get_stop_frequency(), len(data))
        return np.vstack((fpts, data))
    
# =============================================================================
#     def read_data(self):
#         self.instrument.write(':MMEM:STOR:TRAC:DATA TRACE1, \"C:\\Temp\\trace_temp.csv\"')
#         data = self.instrument.query(':MMEM:DATA? \"C:\\Temp\\trace_temp.csv\"')
#         return data
# =============================================================================

    #### Meta
    def take_one(self):
        #print "Acquiring single trace"
        #time.sleep(na.query_sleep*2)
        #na.set_format()
        #time.sleep(na.query_sleep)
        self.clear_averages()
        self.trigger_single()
        time.sleep(self.get_query_sleep())
        self.wait_for_completion()
        ans = self.read_data()
        return ans

    def get_settings(self):
        settings = {"start": self.get_start_frequency(), "stop": self.get_stop_frequency(),
                    "power": self.get_power(), "ifbw": self.get_ifbw(),
                    "sweep_points": self.get_sweep_points(),
                    "averaging": self.get_average_state(), "averages": self.get_averages()
        }
        return settings

    def configure(self, start=None, stop=None, center=None, span=None, resbw=None, vidbw=None, sweep_pts=None,
                  avgs=None, defaults=False, remote=False):
        if defaults:       self.set_default_state()
        if remote:                          self.set_remote_state()
        if start is not None:            self.set_start_frequency(start)
        if stop is not None:            self.set_stop_frequency(stop)
        if center is not None:          self.set_center_frequency(center)
        if span is not None:            self.set_span(span)
        if resbw is not None:         self.set_resbw(resbw)
        if vidbw is not None:            self.set_vidbw(vidbw)
        if sweep_pts is not None:   self.set_sweep_points(sweep_pts)
        if avgs is not None:            self.set_averages(avgs)

    def set_remote_state(self):
        self.set_trigger_source('EXT')
        self.set_timeout(10)

    def set_default_state(self):
        self.set_sweep_points()
        self.set_trigger_source()
        self.instrument.write(":INIT:CONT ON")

    def play_sound(self, freq=None, music=False):
        import winsound
        import time
        soundfile = "tone_files/1kHz_44100Hz_16bit_30sec.wav"
        winsound.PlaySound(soundfile, winsound.SND_FILENAME|winsound.SND_ASYNC)
        time.sleep(1.5)
        # play the system exit sound if set
        winsound.PlaySound("SystemExit", winsound.SND_ALIAS)
        
# =============================================================================
# # In[2]: test instrument
# 
# if __name__ == '__main__':
#     #    condense_nwa_files(r'C:\\Users\\dave\\Documents\\My Dropbox\\UofC\\code\\004 - test temperature sweep\\sweep data','C:\\Users\\dave\\Documents\\My Dropbox\\UofC\\code\\004 - test temperature sweep\\sweep data\\test')
#     sa = E4440("E4440", address="192.168.14.152")
#     print(sa.get_id())
#     print("Taking data")
#     data = sa.take_one()
#     sa.set_default_state()
#     #print data
#     from matplotlib.pyplot import *
#     #data=sa.read_data()
#     plot(data[0], data[1])
#     show()
#     #print "Setting window"
# 
#     #from guiqwt.pyplot import *
#     #nwa_test2(na)
# #    nwa_test2(na)
# #nwa_test2(na)
# 
# # In[]
# import matplotlib.pyplot as plt
# 
# sa = E4440()
# 
# print(sa.get_id())
# 
# data = sa.read_data()
# plt.figure()
# plt.plot(data[0], data[1])
# plt.show()
# 
# # In[]
# import pyvisa
# import matplotlib.pyplot as plt
# 
# RM =  pyvisa.ResourceManager()
# sa = RM.open_resource('TCPIP::192.168.0.60::INSTR')   
# 
# data_str = sa.query(":TRAC:DATA? TRACE1")
# data_conv= np.fromstring(data_str, dtype=float, sep=',')
# 
# plt.figure()
# plt.plot(data_conv)
# plt.show()
# # 
# 
# #    nwa_test3(na)
# 
# 
# =============================================================================
