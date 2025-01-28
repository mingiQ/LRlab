# In[0] : necessary packages
import sys
import  numpy as np
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import SerialInstrument



## DAC LR 
class DAC_CF(SerialInstrument):
    '''
    !!! different from BPF, PySerial requires '\n' at the end of the command!!!!! (2023/5/9)
        /*    ------ communication commands ---------
     *     serial @ 115200, line end character: \n
     *    "?" or "*IDN?"             ->  devstamp         - equipment ID  
     *    "DD dac"                   -> "DD"              - initialize dac
     *    "WF dac val val"           -> "WF"              - set single dac value (float)
     *    "WD dac dval"              -> "WD"              - set single dac value (digital), always in dc mode
     *    "WH dac val val"           -> "WH"              - set double dac & dac+1 to val
     *    "RF dac"                   -> "RF val"          - read dac value (float)
     *    "RD dac"                   -> "RD dval"         - read dac value (digital)
     *    "RH dac"                   -> "RH val"          - read double-dac value
     *    
     *    "WS dac val"               -> "WS"              - write dc change rate [V/s], 0=no ramp
     *    "RS dac"                   -> "RS val"          - read  dc change rate [V/s], 0=no ramp
     *    
     *    "WA dac ac                 -> "WA ac"           - set amplitude of ac modulation in V
     *    "WT dac freq               -> "WT"              - set freq of ac modulation in [Hz]
     *    "WP dac phase              -> "WP"              - set phase of ac modulation in [deg]
     *    "WU dac duty               -> "WU"              - set pulse duty cycle (0-1)
     *    "WM dac mode               -> "WM"              - set modulation type 0-off, 1-sin, 2-square, 3-pulse
     *    "RA dac                    -> "RA float"        - read amplitude of ac modulation in Vrms
     *    "RT dac                    -> "RT float"        - read period of ac modulation in [Hz], 0-disabled
     *    "RP dac                    -> "RP float"        - read phase of ac modulation in [deg]
     *    "RU dac                    -> "RU float"        - read pulse duty cycle
     *    "RM dac                    -> "RM int"          - read modulation type 0-off, 1-sin, 2-square, 3-sin(fast)
     *    "RW dac                    -> "RW many"         - read modulation params "dc,ac,f,ph,du'
     *    
     *    "SR dac"                   -> "SR bits neg pos step dc ldac" - read dac info from ROM
     *    "SW dac neg pos bits"      -> "SW"              - write dac info to ROM
     *    "ER dac"                   -> "ER neg pos bits" - read dac info from EEPROM
     *    "SE dac"                   -> "SE"              - write data from ROM to EEPROM
     *    "ES dac"                   -> "ES"              - write data from EEPROM to ROM
     *    "EC #"                     -> "EC"              - TFT write echo level, 0 (off) 1 (+values) 2 (+command) 3 (+debug)
     *                                                          -1 supresses com responce for WF & WD
     *    
     *    "IC dac Vret Vmin Vmax Vthp Vthn t1 t2 t3 tramp repeats AIN" 
     *                               -> "IC Vres"         - set Ic sweep, returns V step size. AIN > 0(ADC) or 0 (comparator)
     *    "RC dac"                   -> "RC Vret Vmin Vmax Vthp Vthn t1 t2 t3 tramp repeats AIN Vres" 
     *    "IS dac mode"              -> "IS stat/stream"  - performs Ic sweeps, mode=0: <Ic>, std, min, max, repeats
     *                                                                               1: stream of Isw's
     *                                                                               2: strem (Vdac,Vadc) for each step between Vmin and Vmax
     *    "IT dac"                   -> "IT time1, time2" - returns max time in sec for n and 1 sweeps (for IS & IR commands)
     *    "SC AIN+, AIN- event"      -> "SC"              - enables  comparator CMP1
     *    "DC"                       -> "DC"              - disables comparator CMP1
     *    "DI"                       -> "DI"              - enable internal DAC0 (12bit, 3.3V)
     *    "DW val"                   -> "DW"              - set DAC0 val
     *    "PW val"                   -> "PW"              - set pin 23 0/1 LOW/HIGH
     *    
     *    "AR AIN"                   -> "AR val"          - reads value on ADC AIN
     *    "AS ave, cspeed, sspeed, res"  -> "AS rate"     - sets teensy ADC parameters, returns rate [us/read]
     *    
     *    "TS dac dval"              -> "TS"              - comm test: write xAAAAA to dac dval times
     *    "TX pin dval"              -> "TX"              - pin test: write HIGH/LOW to pin dval times
     *    "TT pin"                   -> "TT"              - toggle pin
     *    "ZZ"                       -> "ZZ txt"          - testing output 
     *    
     *    "IRQ"                                           - sets IRQready=false at a serial interrupt level
    '''
    def __init__(self, name='LRdac', address='COM4', enabled=True, timeout=0.1, baudrate=115200):
        SerialInstrument.__init__(self, name, address, enabled, timeout, baudrate)
        self.flush_input()
        self.flush_output()
    
    def get_idn(self):
        self.write('*IDN?')
        return self.read(timeout=self.timeout)
    
    def dac_init(self, ch):
        self.write(f"DD {ch}")
    
    def set_voltage(self, ch, volt):
        self.write(f"DD {ch}")
        self.write(f"WF {ch} {volt} {volt}")
        data = self.read(timeout=self.timeout)
        return data
    
    def get_voltage(self, ch):
        self.write(f"RF {ch}")
        voltage = self.read(timeout=self.timeout)
        #return np.float64(voltage.split('\r\n')[0].split('RF ')[1])
        return np.float64(voltage.split('RF')[1])
        
    
    def set_ramping_rate(self, ch, rate):
        self.write(f"WS {ch} {rate}")
    
    def get_ramping_rate(self, ch):
        self.write(f"RS {ch}")
        rate = self.read(timeout=self.timeout)
        return np.float64(rate.split('\r\n')[0].split('RS ')[1])
        
    def set_ac_mod(self, ch, amp, freq, phase, duty_cycle, mode):
        '''
      "WA dac ac                 -> "WA ac"           - set amplitude of ac modulation in V
     *    "WT dac freq               -> "WT"              - set freq of ac modulation in [Hz]
     *    "WP dac phase 
     *    "WU dac duty               -> "WU"              - set pulse duty cycle (0-1)
     *    "WM dac mode               -> "WM"              - set modulation type 0-off, 1-sin, 2-square, 3-pulse
        '''
        self.write(f"WA {ch} {amp}")
        self.write(f"WT {ch} {freq}")
        self.write(f"WP {ch} {phase}")
        self.write(f"WU {ch} {duty_cycle}")
        self.write(f"WM {ch} {mode}")
        
    def set_duty_cycle(self, ch, duty_cycle):
        self.write(f"WU {ch} {duty_cycle}")
    
    def set_mod_type(self, ch, modtype):
        self.write(f"WM {ch} {modtype}")
        
    def get_ac_mod(self, ch):
        '''
      *    "RW dac                    -> "RW many"         - read modulation params "dc,ac,f,ph,du'
        '''
        config = self.write(f"RW {ch}")
        dc = config.split('DC')[1].split(' AC')[0]
        ac = config.split('DC')[1].split(' AC')[1].split(' F')[0]
        freq = config.split('DC')[1].split(' AC')[1].split(' F')[1].split(' P')[0]
        phase = config.split('DC')[1].split(' AC')[1].split(' F')[1].split(' P')[1].split(' D')[0]
        duty_cycle = config.split('DC')[1].split(' AC')[1].split(' F')[1].split(' P')[1].split(' D')[1].split(' M')[0]
        mod_type = config.split('DC')[1].split(' AC')[1].split(' F')[1].split(' P')[1].split(' D')[1].split(' M')[1]
        
        return np.float64([dc, ac, freq, phase, duty_cycle, mod_type])
    
    
    ###
    ### Experiments
    ###
    def set_IcSweeps(self, ch, vret, vmin, vmax, vthp, vthn, t1, t2, t3, tramp, repeats, adc_ch ):
        '''
        IC dac Vret Vmin Vmax Vthp Vthn t1 t2 t3 tramp repeats AIN
        
        high speed IV sweep for Ic characterization

        Parameters
        ----------
         : TYPE
            dac: dac channel
            Vret: voltage to rest
            Vmin: minimum voltage for sweep
            Vmax: max voltage for sweep
            Vthp: V_threshold, positive (slightly higher than voltage converted + critical current )
            Vthn: V_threshold, negative (slightly lower than voltage converted - critical current)
            t1: stablizing time for DAC
            t2: wait time for DAC
            t3: measurement delay for ADC
            tramp: ramping time
            repeats: number of experiment
            AIN: ADC channel, 0: comparator

        Returns
        -------
        None.

        '''
        self.write(f"IC {ch} {vret} {vmin} {vmax} {vthp} {vthn} {t1} {t2} {t3} {tramp} {repeats} {adc_ch}")
        
    def get_IcSweep_params(self, ch):
        '''
        

        Parameters
        ----------
        ch : assigned dac channel

        Returns
        -------
        decode_param : np.array
        Vret Vmin Vmax Vthp Vthn t1 t2 t3 tramp repeats AIN Vres

        '''
        self.write(f"RC {ch}")
        params = self.read(timeout=self.timeout)
        decode_param = np.float64(params.split('RC')[1].split(' ')[1:])
        return decode_param
        
    def Ic_Sweep(self, ch, sweep_mode):
        '''
        "IS dac mode"              
        -> "IS stat/stream"  
        - performs Ic sweeps, mode=0: <Ic>, std, min, max, repeats
                                   1: stream of Isw's
                                   2: strem (Vdac,Vadc) for each step between Vmin and Vmax
        '''
        self.write(f"IS {ch} {sweep_mode}")
    
# In[]:
    
#dac = DAC_CF()


# In[test aidan dac @ heliox]

#dac = DAC_CF(name='dac', address='COM7', enabled=True, timeout=0.1, baudrate=115200)

# In[]

#dac.__del__()
    
# In[]
    
# =============================================================================
#     
#     def __init__(self, addr):       # COM4
#         self.DAC = serial.Serial(port=str(addr), baudrate=115200, timeout=None)
#         self.DAC.flushInput()
#         self.DAC.flushOutput()
#     def write_read(self, x):                    # CAUTION: all the command should end with '\n'
#         self.DAC.write(bytes( x +'\n', 'utf-8'))
#         time.sleep(0.05)
#         #data = self.arduino.readline().decode('utf-8').rstrip()
#         data = self.DAC.readline().rstrip()
#         print(data)
#     def idn(self):
#         self.DAC.write(bytes('*IDN?\n', 'utf-8'))
#         time.sleep(0.05)
#         data = self.DAC.readline().rstrip()
#         print(data)
#     def voltsweep(self, dac_CF, volt):
#         self.DAC.write(bytes("WF "+ str(dac_CF) +" "+ str(volt) +" "+ str(volt)+"\n", 'utf-8'))
#         time.sleep(0.05)
#         data = self.DAC.readline().rstrip()
#         print(data)
#     def quit_dac(self):
#         self.DAC.close()
# =============================================================================
