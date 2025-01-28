# In[0] : necessary packages
import sys
import numpy as np
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument

import time
import numpy as np
import glob
import os.path

## Copper Mountain VNA


def polar2mag(xs, ys):
    return 20*np.log10(np.sqrt(xs ** 2 + ys ** 2)), np.arctan2(ys, xs) * 180/np.pi


class CMT(VisaInstrument):
    def __init__(self, name='CopperMountain', address='TCPIP0::127.0.0.1::5025::SOCKET', enabled=True, timeout=0.1, rm='C:/Windows/System32/visa64.dll'):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=1
        self.recv_length=65536
        self.term_char=''
        self.instrument.read_termination='\n'
        self.instrument.timeout = 10000
        
    def get_idn(self):
        self.instrument.timeput = 15000
        return self.instrument.query('*IDN?')
     
    def get_query_sleep(self):
        return self.query_sleep
    
     #### frequency setup
    
    def set_start_frequency(self, p):
        self.instrument.write(':SENS:FREQ:STAR ' + str(p))
    
        
    def get_start_frequency(self):
        self.instrument.timeput = 15000000
        return int(np.float64(self.instrument.query(':SENS:FREQ:STAR? ')))

        
    def set_stop_frequency(self, p):
        self.instrument.write(':SENS:FREQ:STOP ' + str(p))
        
        
    def get_stop_frequency(self):
        return int(np.float64(self.instrument.query(':SENS:FREQ:STOP? ')))
     
        
    def set_center_frequency(self, p):
        self.instrument.write(':SENS:FREQ:CENT ' + str(p))
        self.instrument.baud_rate = 57600
         
    def get_center_frequency(self):
        return int(np.float64(self.instrument.query(':SENS:FREQ:CENT? ')))
 
        
    def set_span(self, span):
        self.instrument.write(':SENS:FREQ:SPAN %f' % span)
        self.instrument.baud_rate = 57600
        
    def get_span(self):
        return int(np.int64(self.instrument.query(':SENS:FREQ:SPAN? ')))
 
        
    def set_sweep_points(self, numpts=1600):  
        self.instrument.timeout=160000
        self.instrument.write(':SENSe:SWEep:POINts %d' % numpts)
        
        
    def get_sweep_points(self):
        self.instrument.timeout=160000
        return int(np.float64(self.instrument.query(':SENSe:SWEep:POINts? ')))
    
    
    def set_sweep_mode(self, mode='CONT'):
        """
        Sets the number of trigger signals the specified channel will ACCEPT

        HOLD - channel will not trigger
        CONTinuous - channel triggers indefinitely
        GROups - channel accepts the number of triggers specified with the last SENS:SWE:GRO:COUN <num>.
        SINGle - channel accepts ONE trigger, then goes to HOLD.
        """
        
        self.instrument.write(':INIT ')
        self.instrument.write(':TRIG:SOUR BUS ')

        allowed_modes = ['cont', 'continuous', 'hold', 'sing', 'single']
        if mode.lower() not in allowed_modes:
            return
        elif mode.lower() == 'single':
            self.instrument.write(":TRIG:SEQ:SING ")
        elif mode.lower() == 'hold':
            self.instrument.write(':INIT:CONT OFF')
        elif mode.lower() == 'cont':
            self.instrument.write(':INIT:CONT ON')
            
##    
      
    #### Averaging
    def set_averages(self, averages, channel=1):
        self.instrument.write(":SENS%d:AVERage:COUNt %d" % (channel, averages))

    def get_averages(self, channel=1):
        return int(self.instrument.query(":SENS%d:average:count?" % channel))

    def set_average_state(self, state=True, channel=1):
        if state:
            s = "ON"
        else:
            s = "OFF"
        self.instrument.write(":SENS%d:AVERage:state %s" % (channel, s))

    def get_average_state(self, channel=1):
        return bool(self.instrument.query(":SENS%d:average:state?" % channel))

    def clear_averages(self, channel=1):
        self.instrument.write(":SENS%d:average:clear" % channel)

    def set_ifbw(self, bw, channel=1):
        self.instrument.write("sens%d:bwid %f" % (channel, bw))

    def get_ifbw(self, channel=1):
        return float(self.instrument.query("SENS%d:bwid?" % (channel)))

    def get_operation_completion(self):
        data = self.instrument.query("*OPC?")
        if data is None:
            return False
        else:
            return bool(int(data.strip()))
        
    def set_average_trig(self, trig=True):
        if trig:
            t = 1
        else:
            t = 0
        self.instrument.write(f"TRIG:AVER {t}")
        
    def get_average_trig(self):
        return bool(self.instrument.query(":TRIG:AVER?"))

    def set_trigger_continuous(self, state=True):
        """
        This command sets the trigger mode to continuous (internal) or manual
        NB: to refresh the display, use set_sweep_mode("CONT") in combination
        with this command.
        same as set sweep mode (in CMT)
        """
        if state:
            _state = "on"
        else:
            _state = "off"
        self.instrument.write('initiate:continuous ' + _state)

        
     #### Source

    def set_power(self, power, channel=1):
        self.instrument.write(":SOUR:POW %20.15e dBm " % power)

    def get_power(self, channel=1):
     
        return float(self.instrument.query(":SOUR%d:POW? " % channel))

    def set_output(self, state=True):
        if state or str(state).upper() == 'ON':
            self.instrument.write(":OUTPUT ON")
        elif state == False or str(state).upper() == 'OFF':
            self.instrument.write(":OUTPUT OFF")

    def get_output(self):
        return bool(self.instrument.query(":OUTPUT?"))


    def define_measurement(self, mode, trace=1, channel=1):
        self.instrument.write('CALC%d:PAR%d:DEF %s' % (channel, trace, mode))

    def get_measurements(self, trace=1, channel=1):
        data = self.instrument.query('CALC%d:PAR%d:DEF? ' % (channel, trace))
        return data

    def select_measurement(self, trace=1, channel=1):
        self.instrument.write("calc%d:par%d:sel " % (channel, trace))

    def auto_scale(self, trace=None):
        """
        Performs an Autoscale on the specified trace in the specified window, providing the best fit display.
        Autoscale is performed only when the command is sent; it does NOT keep the trace autoscaled indefinitely.
        """
        if trace is None:
            query = "DISP:WIND:TRAC:Y:AUTO"
        else:
            query = "disp:wind:trac%d:Y:AUTO" % trace
        self.instrument.write(query)
        
        
    #### Trace operation
    
    def set_data_transfer_format(self, format='ascii'):
        """
        Sets the data format for transferring measurement data and frequency data.
        See the Format Commands help section for more help on this topic.
        :param format: Either 'ascii' or 'binary'
        :return:
        """
        send_data = 'ASC,0' if format.lower() == 'ascii' else 'REAL32'
        self.instrument.write("FORM:DATA %s" % send_data)
        if send_data == 'REAL32':
            self._set_byte_order('SWAP')

    def get_data_transfer_format(self):
        """
        Returns the data format for transferring measurement data and frequency data.
        :return: 'ascii' or 'binary'
        """
        self.instrument.timeout = 1600000
        answer = self.instrument.query('FORM:DATA?')
        ret = 'ascii' if 'ASC' in answer else 'binary'
        return ret
    
    def _set_byte_order(self, order='SWAP'):
        """
        #NOTE for Plutonium, the byte order needs to be swapped!
        Set the byte order used for GPIB data transfer. Some computers read data from the analyzer in the reverse order.
        This command is only implemented if FORMAT:DATA is set to :REAL. If FORMAT:DATA is set to :ASCII, the swapped command is ignored.
        :param order: 'swap' for swapped or 'norm' normal order
        :return: None
        """
        if order.upper() in ['SWAP', 'NORM']:
            self.instrument.write("FORM:BORD %s" % order.upper())

    def _get_byte_order(self):
        """
        Returns the byte order used for GPIB data transfer.
        :return: 'SWAP' (swapped) or 'NORM' (normal order)
        """
        return self.instrument.query("FORM:BORD?").strip()
    
    def set_format(self, trace_format='LMAG', trace=1):
        """
        Sets the display format for the measurement.
        This needs to be run after the active trace is set. The following options are available:
            {MLOG|PHAS|GDEL|SLIN|SLOG|SCOM|SMIT|SADM|PLIN|PLOG|POL|MLIN|SWR|REAL|IMAG|UPH}
        """
        allowed = ['MLOG','PHAS','GDEL','SLIN','SLOG','SCOM','SMIT','SADM','PLIN','PLOG','POL','MLIN','SWR','REAL','IMAG','UPH']
        if trace_format.upper() in allowed:
            self.instrument.write("CALC%d:FORM %s" % (trace, trace_format.upper()))
        else:
            raise ValueError("Specified trace format not allowed. Use %s" % allowed)

    def get_format(self, trace=1):
        """set_format: need to run after active trace is set.
        valid options are
         {MLOG|PHAS|GDEL|SLIN|SLOG|SCOM|SMIT|SADM|PLIN|PLOG|POL|MLIN|SWR|REAL|IMAG|UPH}
        """
        data = self.instrument.query("CALC%d:FORM?" % (trace))
        if data is None:
            return data
        else:
            return data.strip()
    
    def set_electrical_delay(self, seconds, channel=1):
        """
        Sets the electrical delay in seconds
        :param seconds: Electrical delay in seconds
        :param channel: Measurement channel
        :return: None
        """
        query = "calc%d:corr:edel:time %e" % (channel, seconds)
        self.write(query)

    def get_electrical_delay(self, channel=1):
        """
        Returns the electrical delay in seconds
        :param channel: Measurement channel
        :return: Electrical delay in seconds
        """
        query = "calc%d:corr:edel:time?" % channel
        data = self.query(query)
        if data is None:
            return None
        else:
            return float(data.strip())

    def set_phase_offset(self, degrees, channel=1):
        """
        Sets the phase offset for the selected measurement
        CALCulate<Ch>[:SELected]:CORRection:OFFSet:PHASe
        :param degrees: Phase offset in degrees. Choose any number between -360 and 360.
        :param channel: Measurement channel
        :return:
        """
        query = "CALC%d[:SEL]:CORR:OFFS:PHAS %.3f" % (channel, degrees)
        self.write(query)

    def get_phase_offset(self, channel=1):
        """
        Returns the phase offset for the selected measurement
        :param channel: Measurement channel
        :return: Numeric, returned value always in degrees
        """
        query = "CALC%d[:SEL]:CORR:OFFS:PHAS?" % channel
        data = self.query(query)
        if data is None:
            return None
        else:
            return float(data.strip())   
        
    #### File Operations
    def save_file(self, fname):
        self.instrument.write('MMEMORY:STORE:FDATA \"' + fname + '\"')

    def read_data(self, sweep_points=None, channel=1, trace=1, timeout=None, data_format=None):
        """
        Read current NWA Data that is displayed on the screen. Returns TWO numbers per data point for Polar ('POL')
        and Smith Chart ('SMIT') format, see set_format.
        :param sweep_points: number of sweep points (optional, saves time)
        :param channel: measurement channel (optional, saves time)
        :param timeout: timeout in seconds (optional)
        :param data_format: 'binary' or 'ascii' (optional, saves time). If specificied, this must be equal to get_data_transfer_format()
        :return: 2 or 3 column data containing the frequency and data (1 or 2 column).
        """
        if data_format is None or not(data_format in ['binary', 'ascii']):
            data_format = self.get_data_transfer_format()

        if sweep_points is None:
            sweep_points = self.get_sweep_points()

        if timeout is None:
            timeout = self.timeout
            
        self.get_operation_completion()
        
        self.instrument.timeout = 160000  
        fdat = self.instrument.query('CALCulate{}:TRACe{}:DATA:FDATa?'.format(channel, trace)) 
        
        data_str = fdat

        if data_format == 'binary':
            len_data_dig = np.int(data_str[1:2])
            len_data_expected = int(data_str[2: 2+len_data_dig])
            len_data_actual = len(data_str[2 + len_data_dig:-1])
            # It may happen that only part of the message is received. We know that this is the case by checking
            # the checksum. If the received data is too short, just read out again.
            while len_data_actual != len_data_expected:
                data_str += b''.join(self.read_lineb(timeout=timeout))
                len_data_actual = len(data_str[2 + len_data_dig:-1])

        data = np.fromstring(data_str, dtype=float, sep=',') if data_format=='ascii' else np.fromstring(data_str[2+len_data_dig:-1], dtype=np.float32)
        fpts = np.linspace(self.get_start_frequency(), self.get_stop_frequency(), sweep_points)
        if len(data) == 2 * sweep_points:
            data = data.reshape((-1, 2))
            data = data.transpose()
            return np.vstack((fpts, data))
        else:
            return np.vstack((fpts, data))

    #### Meta
    def take(self, sweep_points=None, data_format=None):
        """
        Important:
            the PNA-X needs to be in the following mode
                trigger source:IMMediate,
                format:POLar,
                trigger:CONTinuous ON
        :param sweep_points:
            by taking in a sweep_points parameter, we do not need to query the PNA-X for this
            parameter. This way we can save some time.
        :return:
            either fpts, xs, ys,
            or     fpts, mags.
        """
        self.clear_averages()
        # this is the command that triggers the averaging. Execute right before read data.
        self.set_sweep_points(sweep_points)
        self.set_average_trig()
        self.set_sweep_mode('single')
        data = self.read_data(sweep_points=sweep_points, data_format=data_format)
        return data

# =============================================================================
#     def take_in_mag_phase(self, sweep_points=None, data_format=None):
#         _format = self.get_format()
#          
#         fpts, xs, ys = self.take(sweep_points=sweep_points, data_format=data_format)
#         mags, phases = polar2mag(xs, ys)
#         
#         self.set_format(_format)
#         return fpts, mags, phases
# =============================================================================

    def take_one_in_mag_phase(self, sweep_points=None, data_format=None):
        """
        Takes one averaged trace and return fpts, magnitudes and phases
        :param sweep_points: Sweep points (optional, saves time)
        :param data_format: 'ascii' or 'binary' (optional, saves time)
        :return: fpts, mags, phases
        """
        
        _format = self.get_format()
        self.setup_take()
        fpts, xs, ys = self.take(sweep_points=sweep_points, data_format=data_format)
        mags, phases = polar2mag(xs, ys)
        
        self.set_format(_format)
        return fpts, mags, phases

    def setup_take(self, channel=1, averages=None, averages_state=None):
       
        self.set_format('POL')
        self.set_sweep_mode('single')
        #self.set_trigger_continuous()
        
        if averages is not None:
            self.set_averages(averages, channel)
        elif averages_state is not None:
            self.set_average_state(averages_state, channel)
            self.set_trigger_continuous()

    def setup_measurement(self, trace=1, mode=None):
        if mode is None:
            mode = 'S21'
        self.define_measurement(mode, trace, channel=1)
        #self.display_measurement(trace)
        self.select_measurement(trace)
        #self.set_active_trace(1, 1, True)

    def get_settings(self):
        settings = {"start": self.get_start_frequency(), "stop": self.get_stop_frequency(),
                    "power": self.get_power(), "ifbw": self.get_ifbw(),
                    "sweep_points": self.get_sweep_points(),
                    "averaging": self.get_average_state(), "averages": self.get_averages()
                    }
        return settings

    def configure(self, start=None, stop=None, center=None, span=None,
                  power=None, ifbw=None, sweep_points=None, averages=None):
        if start is not None:      self.set_start_frequency(start)
        if stop is not None:       self.set_stop_frequency(stop)
        if center is not None:     self.set_center_frequency(center)
        if span is not None:       self.set_span(span)
        if power is not None:      self.set_power(power)
        if ifbw is not None:       self.set_ifbw(ifbw)
        if sweep_points is not None:  self.set_sweep_points(sweep_points)
        if averages is not None:    self.set_averages(averages)
        
# =============================================================================
#         
# # In[]
# 
# cmt = CMT()
# print(cmt.get_idn()     ) 
# 
# cmt.set_average_trig(trig=True) 
# # In[]
# 
# print(cmt.get_power())
# 
# # In[]
# cmt.set_averages(averages=20)
# 
# cmt.setup_measurement(mode='S21')
# dat = cmt.take_one_in_mag_phase(sweep_points=10)
# 
# print(dat[1])
# print(dat[2])
# # In[]
# 
# #cmt.set_power(power=-30)
# print(cmt.get_start_frequency())
# print(cmt.get_averages())
# print(cmt.get_average_trig())
# 
# =============================================================================
# =============================================================================
# # In[]
# 
# import pyvisa
# 
# RM =  pyvisa.ResourceManager('C:/Windows/System32/visa32.dll')
# CMT= RM.open_resource('TCPIP0::127.0.0.1::5025::SOCKET')       
# 
# 
# # In[]
# CMT.write(':SOUR:POW -30 dBm \n')
# 
# # In[]
# CMT.write(':SOUR:POW -30 dBm \n')
# 
# =============================================================================
