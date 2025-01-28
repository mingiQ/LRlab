# In[0] : necessary packages
import sys
import time
import numpy as np
import glob
import os.path
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument


## Anritsu MS2038C VNA master  ''' TCPIP[board]::host address[::LAN device name][::INSTR] '''

def polar2mag(xs, ys):
    return 20*np.log10(np.sqrt(xs ** 2 + ys ** 2)), np.arctan2(ys, xs) * 180/np.pi


class MS2038(VisaInstrument):
    
    def __init__(self, name='MS2038c', address='TCPIP::192.168.0.105::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=0.05
        self.recv_length=65536
        self.term_char=''
   
    def get_idn(self):
        return self.instrument.query('*IDN?')
     
    def get_query_sleep(self):
        return self.query_sleep
         
    def set_inst_mode(self, mod): 
        self.instrument.write(':INST[:SEL] ' + str(mod))  # SPA : spectrum analyzer, MWVNA: vector network analyzer
        #self._gpib.baud_rate = 57600
        return(self.instrument.query(':INST?'))
    
    #### frequency setup
    
    def set_start_frequency(self, p):
        self.instrument.timeout=160000
        self.instrument.write(':SENS:FREQ:STAR ' + str(p))

    def get_start_frequency(self):
        return int(self.instrument.query(':SENS:FREQ:STAR? '))

        
    def set_stop_frequency(self, p):
        self.instrument.timeout=160000
        self.instrument.write(':SENS:FREQ:STOP ' + str(p))
        
        
    def get_stop_frequency(self):
        self.instrument.timeout=160000
        return int(self.instrument.query(':SENS:FREQ:STOP? '))
     
        
    def set_center_frequency(self, p):
        self.instrument.timeout=160000
        self.instrument.write(':SENS:FREQ:CENT ' + str(p))
     
         
    def get_center_frequency(self):
        return int(self.instrument.query(':SENS:FREQ:CENT? '))
 
        
    def set_span(self, span):
        self.instrument.timeout=160000
        self.instrument.write(':SENS:FREQ:SPAN %f' % span)
        
    def get_span(self):
        return int(self.instrument.query(':SENS:FREQ:SPAN? '))
 
        
    def set_sweep_points(self, numpts=1600):  
        self.instrument.timeout=160000
        self.instrument.write(':SENSe:SWEep:POINts %d' % numpts)
        
        
    def get_sweep_points(self):
        self.instrument.timeout=160000
        return int(self.instrument.query(':SENSe:SWEep:POINts? '))
    
    
    
    def set_sweep_mode(self, mode):
        self.instrument.write(':INITiate:HOLD OFF')
        
        if mode =='single':
            self.instrument.write(':SENSe:SWEep:TYPE SINGle')
        elif mode == 'cont':
            self.instrument.write(':SENSe:SWEep:TYPE CONTinuous')  # defalut val
        else:
            return

   
    #### Averaging
    
    def set_averages(self, averages):
        self.instrument.write(':SENS:AVER:COUN %d' % averages)
        self.instrument.baud_rate = 57600
        
    def get_averages(self):
        return self.instrument.query(':SENS:AVER:COUN?')
    
    def clear_averages(self):
        self.instrument.write(':SENS:AVER:CLE')
        
    def set_ifbw(self, bw):
        self.instrument.write(':SENS:SWE:IFBW %d' % bw)  # must be one of 100000|50000|20000|10000|5000|2000|1000|500|200|100|50|20|10
    
    def get_ifbw(self):
        return self.instrument.query(':SENS:SWE:IFBW ?')
    
    def get_operation_completion(self):
        data = self.instrument.query("*OPC?")
        self.instrument.timeout = 120000  
        if data is None:
            return False
        else:
            return bool(int(data.strip()))  

    
    def set_trigger_continuous(self, state=True):
        """
        This command sets the trigger mode to continuous (internal) or manual
        NB: to refresh the display, use set_sweep_mode("CONT") in combination
        with this command.
        """
        if state:
            _state = "on"
        else:
            _state = "off"
        self.instrument.write(':INITiate:CONTinuous ' + _state)

    def trigger_single(self, channel=1):
        self.instrument.write(':INITiate[:IMMediate]')

    def set_trigger_source(self, source="immediate"):  # IMMEDIATE, MANUAL, EXTERNAL
        allowed_sources = ['ext', 'imm', 'man', 'immediate', 'external', 'manual']
        if source.lower() not in allowed_sources:
            print("source need to be one of " + ', '.join(allowed_sources))
        self.instrument.write('TRIG:SEQ:SOUR ' + source)

    def get_trigger_source(self):  # INTERNAL, MANUAL, EXTERNAL,BUS
        return self.instrument.query(':TRIG:SEQ:SOUR?').strip()


    
    def get_minimum(self):
        mini = self.instrument.write(':CALCulate:MARKer 1 :MINimum')
        print(mini)
        
        
    ### Source
    
    def set_power(self, power):
        '''
        Description: Sets the power levels.
        Syntax: :SOURce:POWer LOW|HIGH
        :SOURce:POWer?
        Cmd Parameter: <char> [LOW|HIGH]
        Query Response: <char> [LOW|HIGH]
        Range: HIGH: 3 dBm to –3 dBm
        LOW: –15 dBm to –25 dBm
        
        '''    
        if power == 'low':
            po = 'LOW'
        elif power == 'high':
            po = 'HIGH'
            
        self.instrument.write(':SOURce:POWer ' + po)
        
    def get_power(self):
        return self.instrument.query(':SOURce:POWer?')
    
    
    #### Trace operation
    
    def set_data_transfer_format(self, format='ascii'):
        """
        Sets the data format for transferring measurement data and frequency data.
        See the Format Commands help section for more help on this topic.
        :param format: Either 'ascii' or 'binary'
        :return:
        """
        send_data = 'ASC,0' if format.lower() == 'ascii' else 'REAL,32'
        self.instrument.write("FORM %s" % send_data)
        if send_data == 'REAL,32':
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
    
    def set_format(self, trace_format='LMAG', trace=1):
        """
        Sets the display format for the measurement.
        This needs to be run after the active trace is set. The following options are available:
        LMAGnitude|SWR|PHASe|REAL|IMAGinary|GDELay|SMITh|
        ISMith|LM/2|LINPolar|LOGPolar|RIMPedance|IIMPedance
        """
        allowed = ['LMAG','SWR','PHAS','REAL','IMAG','GDEL','SMIT','ISM','LM/2','LINP','LOGP','RIMP','IIMP']
        if trace_format.upper() in allowed:
            self.instrument.write("CALC%d:FORM %s" % (trace, trace_format.upper()))
        else:
            raise ValueError("Specified trace format not allowed. Use %s" % allowed)

    def get_format(self, trace=1):
        """set_format: need to run after active trace is set.
        valid options are
        {LMAG|SWR|PHAS|REAL|IMAG|GDEL|SMIT|ISM|LM/2|LINP|LOGP|RIMP|IIMP}
        """
        data = self.instrument.query("CALC%d:FORM?" % (trace))
        if data is None:
            return data
        else:
            return data.strip()
        
    def set_measure(self, mode='S21', channel=1):
        '''
        Description: Defines the S-parameter for the given trace, <Tr>.
        <Tr> is the trace number in the range 1 to 4. If no trace number is
        specified, then the <Tr> parameter defaults to trace number 1. The
        query version of this command returns “S11” if the S-parameter is set
        to S11, “S21” if set to S21, “S12” if set to S12, “S22” if set to S22,
        “SD1D1” if set to SD1D1, “SC1C1” if set to SC1C1, “SC1D1” if set to
        SC1D1, and “SD1C1” if set to SD1C1.
        Note that S-parameter S D1D1 , S C1C1 , S C1D1 , and S D1C1 are available
        only if option 77 is installed.
        Syntax: [:SENSe]:TRACe<Tr>:SPARams
        S11|S21|S12|S22|SD1D1|SC1C1|SC1D1|SD1C1
        [:SENSe]:TRACe<Tr>:SPARams?
        Cmd Parameter: <char> [S11|S21|S12|S22|SD1D1|SC1C1|SC1D1|SD1C1]
        Query Response: <char> [S11|S21|S12|S22|SD1D1|SC1C1|SC1D1|SD1C1]
        '''
        self.instrument.write('[:SENSe]:TRACe{}:SPARams {}'.format(channel, mode))
        measurement_mode = self.instrument.query('[:SENSe]:TRACe<Tr>:SPARams?')
        print('trace {} is set to measure {}'.format(channel, measurement_mode))
        
    #### File operations
        
    def save_s2p(self, filename): #save s2p file
        self.instrument.write(':MMEMory:MSIS INTernal')                  # VERY important to define temp data storage first!
        self.instrument.timeout = 120000                                 # Anritsu is really lazy.. give it sufficient time to respond...:( 
        self.instrument.write(':MMEM:STOR:TRAC 4, \"tempdata\"')
        self.instrument.baud_rate = 57600
        time.sleep(1)
        self.instrument.timeout = 120000
        data = self.instrument.query(':MMEM:DATA? \"tempdata.s2p\"')
        self.instrument.baud_rate = 57600
        self.instrument.timeout = 120000
        time.sleep(2)
        data = data.replace('\r','')
        data = data[data.index('!'):]
        time.sleep(1)
        with open(filename,'w') as f:
            f.write(data[data.index('!'):])
            
    def read_data(self, sweep_points=None, channel=1, timeout=None, data_format=None):
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
        #self.get_operation_completion()
        #self.read(timeout=0.1)
        self.instrument.timeout = 160000             # anritsu is lazy - apply sufficient timeout!
        fdat = self.instrument.query("CALC%d:DATA? FDAT" % channel)
        #time.sleep(0.1*sweep_points)
        
        data_str = fdat[8:]

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
            the Anritsu needs to be in the following mode
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
        self.set_sweep_mode('single')
        time.sleep(0.01*sweep_points)
        data = self.read_data(sweep_points=sweep_points, data_format=data_format)
        return data
    
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
        mags, phases = xs, ys
       
        self.set_format(_format)
        return fpts, mags, phases
    
    def setup_take(self, averages=None):
      
        self.set_format('LOGP')
        self.set_sweep_mode('single')

        if averages is not None:
            self.set_averages(averages, True)
    


    def clear_traces(self):
        self.delete_trace()
        self.delete_measurement()

    def setup_measurement(self, name, mode=None):
        if mode is None:
            mode = name
        self.define_measurement(name, 1, mode)
        self.display_measurement(name)
        self.select_measurement(name)
        self.set_active_trace(1, 1, True)

    def get_settings(self):
        settings = {"start": self.get_start_frequency(), "stop": self.get_stop_frequency(),
                    "power": self.get_power(), "ifbw": self.get_ifbw(),
                    "sweep_points": self.get_sweep_points(),
                     "averages": self.get_averages()
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
        if averages is not None:      self.set_averages(averages)
            
    def save_trace(self,filename):
        temp = time.strftime("%d%b%Y%H%M%S", time.localtime())
        self.instrument.write(':MMEMory:MSIS INTernal')
        self.instrument.timeout = 120000     
        self.instrument.write(':MMEM:STOR:TRAC 0,\"'+str(temp)+'\"')
        self.instrument.baud_rate = 57600
        self.instrument.timeout = 120000   
        data = self.instrument.query(':MMEMory:DATA? \"'+str(temp)+'.spa\"')
        self.instrument.baud_rate = 57600
        self.instrument.timeout = 120000
        data = data.replace('\r','')
        data = data.replace('MHz', '')
        data = data.replace('=', ',')
        data = data[data.index('\nP_'):]
        data = data[:data.index('#')]
        with open(filename,'w') as f:
            f.write(data.replace('P_', ''))
            
    def transfer(self, filename, path): # transfer internal memory to pc
        data = self.instrument.query(':MMEM:DATA? \"'+ filename +'\"')
        data = data.replace('\r','')
        data = data[data.index('!'):]
        with open(path+filename,'w') as f:
            f.write(data[data.index('!'):])
    
# =============================================================================
# # In[]
# 
# vna = MS2038(address='TCPIP::192.168.0.110::INSTR')   
# print(vna.get_idn())    
# 
# # In[]
# 
# a = vna.get_data_transfer_format()
# print(a)
# 
# a = vna.read_data(timeout=160000)   
# 
# # In[]
# a = vna.take_one_in_mag_phase(sweep_points=1000)    
# 
# # In[]
# 
# vna.configure(start=None , stop=None, center=4e9, span=60e6,
#                   power='low', ifbw=1e3, sweep_points=2001, averages=5)
# 
# #vna.setup_measurement(mode="S21")
# =============================================================================
