# In[0] : necessary packages
import sys
import numpy as np
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument

import time
import numpy as np
import glob
import os.path

## Hantek Arbitrary waveform generator


def polar2mag(xs, ys):
    return 20*np.log10(np.sqrt(xs ** 2 + ys ** 2)), np.arctan2(ys, xs) * 180/np.pi

'''
USB0::0x049F::0x505B::HTG10000522222::INSTR
'''

class HDG2062(VisaInstrument):
    def __init__(self, name='HantekAWG', address='USB0::0x049F::0x505B::HTG10000522222::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=1
        self.recv_length=65536
        self.term_char=''
        #self.instrument.read_termination='\n'
        self.instrument.timeout = 160000
        
    def get_idn(self):
        #self.instrument.timeout = 15000
        self.instrument.write('*CLS')
        return self.instrument.query('*IDN?')
    
    def get_operation_completion(self):
        data = self.instrument.query("*OPC?")
        if data is None:
            return False
        else:
            return bool(int(data.strip()))
        
    def set_output(self, ch, state=True):
        """
        Disable or enable the Output connector on the front panel. The default is ON. The Output key is lit when enabled.
        :param state: bool
        :return: None
        """
        if state:
            self.instrument.timeout = 160000
            self.instrument.write('*CLS')
            self.instrument.write(f'OUTPUT{ch} ON')
            self.instrument.write('*CLS')
        else:
            self.instrument.timeout = 160000
            self.instrument.write('*CLS')
            self.instrument.write(f'OUTPUT{ch} OFF')
            self.instrument.write('*CLS')

    def get_output(self, ch):
        """
        “0” or “1” indicating the on/off state of the Output connector on the front panel is returned.
        :return: bool
        """
        self.instrument.write('*CLS')
        status = self.instrument.query(f'OUTPUT{ch}?')
        output = {'ON':1, 'OFF':0}
        return output[status]

    def set_output_polarity(self, ch, state='normal'):
        """
        Invert the waveform relative to the offset voltage. The default is NORM, in which
        the waveform goes positive during the first part of the cycle and in INV mode the
        waveform goes negative during the first part of the cycle. The offset remains the
        same when the waveform is inverted and the Sync signal is not inverted.
        :param state: string, 'normal' or 'inverted'
        :return: None
        """
        if state.lower() in ['normal', 'inverted']:
            do_set = 'POS' if state.lower() == 'normal' else 'NEG'
            self.instrument.write(f'SOUR{ch}:MOD:PSK:POL {do_set}')
        else:
            print("State must be 'normal' or 'inverted' for inverted output.")

    def get_output_polarity(self, ch):
        """
        Query the polarity of the waveform. “NORM” or “INV” indicating the polarity will be returned.
        :return: string
        """
        answer = self.instrument.query(f'SOURce{ch}:MOD:PSKey:POLarity?')
        return answer.strip()

# =============================================================================
#     def set_termination(self, load=None):
#         """
#         Select the desired output termination. It can be any value (in ohms) between 1Ω
#         and 10kΩ. INF sets the output termination to “high impedance” (>10 kΩ). The
#         default is 50Ω. The specified value is used for amplitude, offset, and high/low
#         level settings.
#         :param load: float or string
#         :return: None
#         """
#         if load: self.write('OUTPUT:LOAD %s' % load)
# 
#     def get_termination(self):
#         """
#         Query the current load setting in ohms. Return the current load setting or “9.9E+37” meaning “high impedance”.
#         :return: float
#         """
#         return float(self.instrument.query('OUTPUT:LOAD?'))
# =============================================================================

    def set_function(self, ch, ftype="sine"):
        """
        Sets the function. ftype must be one of the strings 
        SINusoid|SQUare|RAMP|
        PULSe|NOISe|DC|SINC|EXPFall|HAVErsine
        |LOREntz| DUALtone|GAUSe|
        ECG| USER| HARMonic|, default SINusoid
        """
        ftypes = {"SINE": "SIN", 
                  "SQUARE": "SQU", 
                  "SQU": "SQU", 
                  "RAMP": "RAMP", 
                  "PULSE": "PULS", 
                  "NOISE": "NOIS",
                  "DC": "DC",
                  "SINC" : "SINC",
                  "EXP" : "EXPF",
                  "HAVERSINE" : "HAVE",
                  "LORENTZ" : "LORE",
                  "GAUSE" : "GAUS",
                  "DUAL" : "DUAL",
                  "HARMOINC" : "HARM",
                  "USER": "USER"}
        ftype_str = ftypes[ftype.upper()]
        self.instrument.write('*CLS')
        self.instrument.write(f'SOURce{ch}:FUNCtion {ftype_str}')

    def set_square(self, ch, dutycycle=50):
        """
        There are limitations to the duty cycle you can set. 
        The pulse duty cycle is defined as:
        Duty Cycle = 100 x Pulse Width ÷ Period
        Pulse width is the time from the 50% threshold of a pulse's rising edge to the 50% threshold of
        the next falling edge.
        The pulse duty cycle range is 0 percent to 100 percent. However, the pulse duty cycle is limited
        by minimum pulse width and edge time restrictions, which prevent you from setting exactly 0
        percent or 100 percent. For example, for a 1 kHz pulse waveform, you are typically restricted to
        pulse duty cycles in the range 0.002 percent to 99.998 percent, limited by the minimum pulse
        width of 16 ns.
        Restrictions Based on Pulse Width: The specified pulse duty cycle must conform to the following
        restrictions determined by the minimum pulse width (Wmin).
        Duty Cycle ≥ (Wmin / Period) X 100
        Duty Cycle ≤ (1 – Wmin / Period) X 100
        """
        self.instrument.write('*CLS')
        self.instrument.write(f'SOURce{ch}:FUNCtion:SQUare:DCYCle {dutycycle}')
        #return self.instrument.write('FUNC:SQU:DCYC %s' % str(dutycycle))

    def get_function(self, ch):
        """
        Query the selection made by FUNCtion USER command. One of the strings “SIN”,
        “SQU”, “RAMP”, “PULS”, “NOIS”, “DC”, and “USER” will be returned.
        :return: string
        """
        self.instrument.write('*CLS')
        func = self.instrument.query(f'SOURce{ch}:FUNCtion?')
        return func

    def set_frequency(self, ch, frequency):
        """
        Set the frequency in Hz
        :param frequency: float
        :return: None
        """
        self.instrument.write('*CLS')
        self.instrument.write(f'SOURce{ch}:FREQuency {frequency}')

    def get_frequency(self, ch):
        """
        Returns the frequency in Hz
        :return: float
        """
        self.instrument.write('*CLS')
        return float(self.instrument.query(f'SOURce{ch}:FREQuency?'))
# =============================================================================
# 
#     def set_period(self, period):
#         """
#         Specify the pulse period. The range is from 200 ns to 2000 seconds. The default is 1 ms.
#         :param period: float
#         :return: None
#         """
#         self.write('PULS:PER %f' % (float(period)))
# 
#     def get_period(self):
#         """
#         The period of the pulse waveform will be returned in seconds.
#         :return: float
#         """
#         return float(self.instrument.query('PULS:PER?'))
# 
# =============================================================================
    def set_pulse_width(self, ch, width):
        """
        Specify the pulse width in seconds. The range is from 16 ns to 1000 us. The default is 500μs.
        :param width: float (seconds)
        :return: None
        """
        self.instrument.write('*CLS')
        self.instrument.write(f'SOURce{ch}:FUNCtion:PULSe:WIDTh {width}')

    def get_pulse_width(self, ch):
        """
        Query the pulse width. The pulse width in seconds will be returned.
        :return: float
        """
        self.instrument.write('*CLS')
        return float(self.instrument.query(f'SOURce{ch}:FUNCtion:PULSe:WIDTh?'))

    def set_pulse_duty_cycle(self, ch, percent):
        """
        The default is 50 percent. 
        
        The pulse duty cycle is defined as:
        Duty Cycle = 100 x Pulse Width ÷ Period
        Pulse width is the time from the 50% threshold of a pulse's rising edge to the 50% threshold of
        the next falling edge.
        The pulse duty cycle range is 0 percent to 100 percent. However, the pulse duty cycle is limited
        by minimum pulse width and edge time restrictions, which prevent you from setting exactly 0
        percent or 100 percent. For example, for a 1 kHz pulse waveform, you are typically restricted to
        pulse duty cycles in the range 0.002 percent to 99.998 percent, limited by the minimum pulse
        width of 16 ns.
        Restrictions Based on Pulse Width: The specified pulse duty cycle must conform to the following
        restrictions determined by the minimum pulse width (Wmin).
        Duty Cycle ≥ (Wmin / Period) X 100
        Duty Cycle ≤ (1 – Wmin / Period) X 100
        
        :param percent: float
        :return: None
        """
        self.instrument.write('*CLS')
        self.instrument.write(f'SOURce{ch}:FUNCtion:PULSe:DCYCle {percent}')

    def get_pulse_duty_cycle(self, ch):
        """
        Query the pulse duty cycle. The duty cycle in percent will be returned.
        :return: float
        """
        self.instrument.write('*CLS')
        return float(self.instrument.query(f'SOURce{ch}:FUNCtion:PULSe:DCYCle?'))
    
    '''
    voltage subsystem
    ?The relationship between offset voltage and output amplitude is shown below.
    |Voffset| < Vmax - Vpp/2
    ?Setting the high and low levels also sets the waveform amplitude and offset. For example, if you
    set the high level to +2 V and the low level to -3 V, the resulting amplitude is 5 Vpp, with a -500 mV offset)
    '''

    def set_offset(self, ch, offset):
        """
        Specify the dc offset voltage. The default is 0 volts. The minimum value is the
        most negative dc offset for the chosen function and amplitude and the maximum
        value is the largest dc offset for the chosen function and amplitude.
        :param offset: float
        :return: None
        """
        self.instrument.write(f"SOUR{ch}:VOLT:OFFSET {offset}")
    
    def get_offset(self, ch):
        """
        Query the dc offset voltage for the current function.
        :return: float
        """
        return float(self.instrument.query("SOUR{ch}:VOLT:OFFSET?"))
    
    def set_amplitude(self, ch, amp):
        """
           The VOLTage subsystem sets parameters related to output voltage
        """
        self.instrument.write(f"SOURce{ch}:VOLTage {amp}")
        
    def get_anplitude(self, ch):
        
        return float(self.instrument.query(f'SOURce{ch}:VOLTage?'))
        
    
# =============================================================================
#     def set_voltage_high(self, high):
#         """
#         Specify the high voltage level. The default high level for all functions is +50 mV.
#         :param high: float
#         :return: None
#         """
#         self.instrument.write('VOLT:HIGH %.5f' % high)
#     
#     def get_voltage_high(self):
#         """
#         Query the high voltage level.
#         :return: float
#         """
#         answer = self.instrument.query('VOLT:HIGH?')
#         return float(answer.strip())
#     
#     def set_voltage_low(self, low):
#         """
#         Specify the low voltage level. The default low level for all functions is -50 mV.
#         :param high: float
#         :return: None
#         """
#         self.instrument.write('VOLT:LOW %.5f' % low)
#     
#     def get_voltage_low(self):
#         """
#         Query the low voltage level.
#         :return: float
#         """
#         answer = self.instrument.query('VOLT:LOW?')
#         return float(answer.strip())
# 
# =============================================================================

    """
    Burst subsystem
    """
    
    def set_burst_cycles(self, ch, cycles=1):
        """
        Specify the number of cycles to be output in each burst (triggered burst mode
        only). The range is from 1 to 50,000 cycles in 1 cycle increments and the default
        is 1 cycle.
        :param cycles: integer
        :return: None
        """
        

        self.instrument.write(f'SOURce{ch}:BURSt:NCYCles {cycles}')

    def get_burst_cycles(self, ch):
        """
        The burst count will be returned. The range is from 1 to 50,000, and 9.9E+37 is
        returned if Infinite is specified.
        :return: integer
        """
        answer = self.instrument.query(f'SOURce{ch}:BURSt:NCYCles?')
        return int(float(answer))
    
    
    def set_burst_state(self, ch, state=True):
        """
        Disable or enable the burst mode.
        :param state: bool
        :return: None
        """
        if state:
            self.instrument.write(f'SOURce{ch}:BURSt ON')
        else:
            self.instrument.write(f'SOURce{ch}:BURSt OFF')

    def get_burst_state(self, ch):
        """
        “0” (OFF) or ”1” (ON) will be returned.
        :return: integer
        """
        return self.instrument.query(f'SOURce{ch}:BURSt?')
    
    def set_burst_mode(self, ch, mode):
        if mode == 'triggered':
            self.instrument.write(f'SOURce{ch}:BURSt:MODE TRIGgered')
        elif mode == 'gated':
            self.instrument.write(f'SOURce{ch}:BURSt:MODE GATed')
        elif mode == 'infinity':
            self.instrument.write(f'SOURce{ch}:BURSt:MODE INFinity')

    def get_burst_mode(self, ch):
        return self.instrument.query(f'SOURce{ch}:BURSt:MODE?').strip()

# =============================================================================
#     def set_burst_period(self, period):
#         """
#         Specify the burst period for bursts with internal (immediate) trigger source. The
#         burst period is ignored when external or manual trigger source is enabled (or
#         when the gated burst mode is chosen).
#         :param period: period in seconds
#         :return: None
#         """
#         self.instrument.write('BURSt:INTernal:PERiod %f' % period)
# 
#     def get_burst_period(self):
#         """
#         The burst period in seconds will be returned.
#         :return: float
#         """
#         return float(self.instrument.query('BURSt:INTernal:PERiod?'))
# 
# =============================================================================
    def set_burst_phase(self, phase):
        """
        Specify the starting phase in degrees or radians according to UNIT:ANGL
        command. The range is from -360 degrees to +360 degrees (or from -2Π to
        +2Π radians) and the default is 0 degree (0 radians).
        """
        self.instrument.write('BURSt:PHASe %d' % phase)

    def get_burst_phase(self):
        """
        The starting phase in degree or radians will be returned.
        :return: float
        """
        return float(self.instrument.query('BURSt:PHASe?'))




    def set_burst_trigger_slope(self, edge):
        edge = edge.upper();
        if edge == 'POS' or edge == 'POSITIVE':
            self.instrument.write('TRIGger:SLOPe %s' % "POSitive")
        elif edge == 'NEG' or edge == 'NEGATIVE':
            self.instrument.write('TRIGger:SLOPe %s' % "NEGative")

    def get_burst_trigger_slope(self):
        return self.instrument.query('TRIGger:SLOPe?')[:-1]
    
    
    '''
    Display control
    '''
    def set_display_saver(self, state=False):
        self.instrument.write('DISPlay:SAVer ON')


#==============================================================================================
    def set_pulse_transition(self, frequency):
        self.instrument.write('FUNC:PULS:TRAN %f' % (float(frequency)))

    def get_pulse_transition(self):
        return float(self.instrument.query('FUNC:PULS:TRAN?'))


    def set_autorange(self, range):
        """
        Disable or enable the voltage auto-ranging. The default is “On” where the
        waveform generator selects an optimal setting for the output amplifier and
        attenuators.
        :param range: 'On' or 'Off'
        :return: None
        """
        self.instrument.write('VOLT:RANGE:AUTO %s' % range.upper())

    def get_autorange(self):
        """
        “0” (off) or “1” (on) indicating the auto-ranging enable state is returned.
        :return: bool
        """
        return int(self.instrument.query('VOLT:RANGE:AUTO?').split('\n')) == 1

  

    def set_trigger_source(self, source="INT"):
        """
        Specify a trigger source for the triggered burst mode only. The waveform
        generator accepts a software (BUS) trigger, an immediate (internal) trigger, or
        a hardware trigger from the rear-panel EXT TRIG connector. The default is IMM.
        :param source:
        :return: None
        """
        trig_types = {'INT': 'IMM', 'INTERNAL': 'IMM', 'EXTERNAL': 'EXT', 'EXT': 'EXT', 'BUS': 'BUS', 'MAN': 'MAN'}
        trig_type_str = trig_types[source.upper()]
        self.instrument.write('TRIG:SOURCE %s' % trig_type_str)

    def get_trigger_source(self):
        """
        Query the trigger source. “IMM” or “BUS” or “EXT” string indicating the trigger
        source will be returned.
        :return: string
        """
        return self.instrument.query('TRIG:SOURCE?').strip()

    def set_trigger_out(self, state):
        """
        Disable or enable the trigger out signal. The default is OFF. When the trigger out
        signal is enabled, a TTL-compatible square waveform with the specified edge is
        output from the Ext Trig connector on the rear panel at the beginning of the
        sweep or burst.
        :param state: bool
        :return: None
        """
        if state:
            self.instrument.write('OutPut:TRIGger %s' % "ON")
        else:
            self.instrument.write('OutPut:TRIGger %s' % "OFF")

    def get_trigger_out(self):
        """
        “0” or “1” indicating the trigger out signal state will be returned.
        :return: bool
        """
        if self.instrument.query('OutPut:TRIGger?') == '1\n':
            return True
        else:
            return False

    def set_trigger_slope(self, edge):
        """
        Specify an edge for the “trigger out” signal.
        :param edge: string ('POS' or 'NEG')
        :return: None
        """
        edge = edge.upper();
        if edge == 'POS' or edge == 'POSITIVE':
            self.instrument.write('OutPut:TRIGger:SLOPe %s' % "POSitive")
        elif edge == 'NEG' or edge == 'NEGATIVE':
            self.instrument.write('OutPut:TRIGger:SLOPe %s' % "NEGative")

    def get_trigger_slope(self):
        """
        “POS” or “NEG” string indicating the edge for the “trigger out” signal will be returned.
        :return: 'POS' or 'NEG'
        """
        return self.instrument.query('OutPut:TRIGger:SLOPe?')[:-1]

    def trigger(self):
        """
        Issue an immediate trigger from the remote interface. This command can
        trigger a sweep or burst with any available trigger source (TRIG:SOUR
        command).
        :return: None
        """
        self.instrument.write('TRIGGER')


    def set_symmetry(self, percent):
        """
        Specify the symmetry percentage for ramp waves. Symmetry represents the
        amount of time per cycle that the ramp wave is rising (supposing the waveform
        polarity is not inverted). The range is from 0% to 100% and the default is 100%.
        :param percent: float
        :return: None
        """
        self.instrument.write('FUNCtion:RAMP:SYMMetry %d' % percent)

    def get_symmetry(self):
        """
        Query the current symmetry setting in percent.
        :return: integer
        """
        return int(self.instrument.query('FUNCtion:RAMP:SYMMetry?'))

    def send_waveform(self, vData, Vpp=2):
        """Rescale and send waveform data to the Tek"""
        # get range and scale to U16
        vI16 = self.scale_waveform_to_I16(vData, Vpp)
        length = len(vI16)
        # create data as string with header
        sLen = '%d' % (2 * length)
        sHead = ':DATA:DAC VOLATILE, #%d%s' % (len(sLen), sLen)
        # write header + data
        self.instrument.write(sHead + vI16.tostring())
        # select volatile waveform
        self.instrument.write(':FUNC:USER VOLATILE')

    def scale_waveform_to_I16(self, vData, dVpp):
        """Scales the waveform and returns data in a string of I16"""
        # clip waveform and store in-place
        np.clip(vData, -dVpp / 2., dVpp / 2., vData)
        vI16 = np.array(2047 * vData / (dVpp / 2.), dtype=np.int16)
        return vI16

    def get_settings(self):
        settings = SocketInstrument.get_settings(self)
        settings['id'] = self.get_id()
        settings['output'] = self.get_output()
        settings['frequency'] = self.get_frequency()
        settings['function'] = self.get_function()
        settings['amplitude'] = self.get_amplitude()
        settings['offset'] = self.get_offset()
        return settings

    def set_exp_trigger(self):

        self.set_output(False)
        self.set_function('PULSE')
        self.set_frequency(5e4)
        self.set_offset(0.)
        self.set_amplitude(1.5)
        self.instrument.write('FUNCtion:PULSe:TRANsition 5e-9')
        self.write('FUNCtion:PULSe:WIDTh 50e-9')
        self.set_output(True)

    
# In[debug] 
    


    
# In[test]

awg = HDG2062()


# In[]

print(awg.get_idn())

# In[]

a = awg.get_operation_completion()
print(a)
# In[]

awg.set_output(ch=2, state=1)

# In[]
awg.set_output(ch=2, state=0)

# In[]

import pyvisa
rm = pyvisa.ResourceManager()

awg = rm.open_resource('USB0::0x049F::0x505B::HTG10000522222::INSTR')

# In[]
awg.timeout = 16000
a = awg.query('*IDN?\n')
print(a)

# In[]

awg.write('*CLS')

# In[]
awg.write('*CLS')
awg.timeout = 160000
a = awg.query('*OPC?')
print(a)

# In[]
awg.write('*CLS')
awg.timeout = 160000
awg.write('OUTPUT1 ON')
awg.write('*CLS')

# In[]
awg.write('*CLS')
awg.timeout = 16000

awg.write('SOURce1:FUNCtion PULS')

# In[]
awg.write('*CLS')
awg.timeout = 16000
a = awg.query(f'SOURce1:FUNCtion?')
print(a)


# In[]
awg.write('*CLS')
awg.timeout = 16000
awg.write('SOURce1:BURSt ON')
#print(a)

awg.write('*CLS')
awg.timeout = 160000
a = awg.query('SOURce1:BURSt?')
print(a)



# In[]
awg.timeout = 160000
a = awg.query('SOUR1:MOD:PSK:SOUR?')