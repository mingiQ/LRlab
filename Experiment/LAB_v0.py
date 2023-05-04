# In[0] : necessary packages
import numpy as np
import pyvisa
import serial
import socket
import json
import time
import sys
sys.path.append('Z:/Mingi Kim/FPGA/qick/qick_lib')
from qick import QickConfig
import Pyro4

# github desktop test

# In[1]: VISA controlled instruments

rm = pyvisa.ResourceManager()
def numtostr(mystr):
    return '%20.15e' % mystr

# In[2]: instruments

### Leiden temperature monitor
class Thermo_leiden:
    def __init__(self):
        self.HOST = '10.164.26.139'
        self.PORT = 8888
        self.soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.soc.connect((self.HOST, self.PORT))
        self.repl = self.soc.recv(1024)
        print(self.repl.decode())
        
    def readT(self, addr):
        self.cmd = str(addr)                #3K	Still	CP	MC TT-1701	3k-pr	Still-Pr/Crnx/PT	CP-50mK-Pr	MC-pr	CMN 089 MC	CMN 060 pr
        self.soc.send(self.cmd.encode())
        repl = self.soc.recv(1024)
        return float(repl)
    
    def Tcal(self, addr):
        self.cmd = 'CAL' + str(addr)
        self.soc.send(self.cmd.encode())
        repl = self.soc.recv(1024)
        return float(repl)
    
    def quit_Tc(self):
        self.cmd = 'quit'
        self.soc.send(self.cmd.encode())
        
    def close_Tc(self):
        self.cmd = 'close'
        self.soc.send(self.cmd.encode())
        

# =============================================================================
# # In[thermo code debug] will add Tc cal code
# 
# HOST = '10.164.26.139'
# PORT = 8888
# soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
# soc.connect((HOST, PORT))
# repl = soc.recv(1024)
# print(repl.decode())
# 
# # In[]
# msg = '9'
# soc.send(msg.encode())
# 
# # In[]
# 
# repl = soc.recv(1024)
# print(repl.decode())
# 
# # In[]
# =============================================================================

### Keithley digital multimeter 2612
class Keithley_2612:
    def __init__(self, addr):
        self._gpib = rm.open_resource('GPIB0::'+str(addr)+'::INSTR')
    def __del__(self):
        self._gpib.close()
    def setvolt(self,v): # numpy float v
        s = np.array2string(v,precision=6)
        self._gpib.write('smua.source.levelv = '+s)
    def volt_on(self):
        self._gpib.write('smua.source.output =smua.OUTPUT_ON')
    def volt_off(self):
        self._gpib.write('smua.source.output =smua.OUTPUT_OFF')
    def meas_i(self, t):
        self._gpib.write('smua.nvbuffer1.clear()')
        self._gpib.write('smua.measure.i(smua.nvbuffer1)')
        time.sleep(t)
        self._gpib.write('printbuffer(1,1, smua.nvbuffer1.readings)')
        x = self._gpib.read()
        v = np.float_(np.fromstring(x,dtype=float,sep=' '))
        return v      
        
### Keythely nanovoltmeter 2182a
class nanovolt:
    def __init__(self, addr):
        self._gpib = rm.open_resource('GPIB0::'+str(addr)+'::INSTR')
    def volt(self):
        self._gpib.write('FETCH?')
        x = self._gpib.read()
        v = np.float_(np.fromstring(x,dtype=float,sep=' '))
        return v
    # def null(self):
    #     self._gpib.write('FETCH?')
    #     x = self._gpib.read()
    #     v = np.float_(np.fromstring(x,dtype=float,sep=' '))
    #     self._gpib.write(':SENS:VOLT:REF 0:STAT 1' +  str(v) + 'ACQ')
    
    
### Keysight 34461A digital voltmeter
class DMM:
    def __init__(self, addr):
        self._gpib = rm.open_resource('GPIB0::'+str(addr)+'::INSTR')
    def __del__(self):
        self._gpib.close()
    def volt(self):
        self._gpib.write('READ?')
        x = self._gpib.read()
        v = np.float_(np.fromstring(x,dtype=float,sep=' '))
        return v
    def null(self, s):
        self._gpib.write('SENS:VOLT:NULL ' + str(s))
        

## MicroLambda bandpass filter --> need to fix!
class MLBPF:
    def __init__(self, addr):
        self.arduino = serial.Serial(port=str(addr), baudrate=115200, timeout=0.1)
    def write_read(self, x):
        self.arduino.write(bytes(x, 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().decode('utf-8').rstrip()
        print(data)
    def idn(self):
        self.arduino.write(bytes('?', 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().decode('utf-8').rstrip()
        print(data)
    def freq(self, f):
        self.arduino.write(bytes('FH' + ' ' + str(f), 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().decode('utf-8').rstrip()
        print(data)
    def quit_bpf(self):
        self.arduino.close()
        

## AMI430 magnet

class AMI430:
    def __init__(self, addr):       # COM7
        self.magnet = serial.Serial(port=str(addr), baudrate=115200, timeout=.1)
        
        
## DAC LR
class DAC_CF:
    def __init__(self, addr):       # COM4
        self.DAC = serial.Serial(port=str(addr), baudrate=115200, timeout=.1)
    def idn(self):
        self.DAC.write(bytes('*IDN?', 'utf-8'))
        data = self.DAC.readline().decode('utf-8').rstrip()
        print(data)
    def voltsweep(self, dac_CF, volt):
        self.DAC.write(bytes("WF "+ str(dac_CF) +" "+ str(volt) +" "+ str(volt), 'utf-8'))
        time.sleep(0.05)
        data = self.arduino.readline().decode('utf-8').rstrip()
        print(data)
# =============================================================================
# # In[debug]:
# DAC = serial.Serial(port="COM4", baudrate=115200, timeout=0)
# # In[]
# DAC.write(bytes('?', 'utf-8'))
# data = DAC.readline().decode('utf-8').rstrip()
# print(data)
# 
# # In[]:
# DAC.write(bytes("WF 8 1 1", 'utf-8'))
# time.sleep(0.05)
# data = DAC.readline().decode('utf-8').rstrip()
# print(data)
# 
# =============================================================================


    
        
#   "WF dac val val"           -> "WF"              - set single dac value (float)
        
        

## Anritsu MS2038C VNA master  ''' TCPIP[board]::host address[::LAN device name][::INSTR] '''
## TCP/IP = 192.168.0.105
## Tested in / RFstation pc / CF room pc /

class MS2038:
    def __init__(self, addr):
        self._gpib = rm.open_resource('TCPIP::'+str(addr)+'::INSTR')
    def mode(self, mod): 
        self._gpib.write(':INST[:SEL] ' + str(mod))  # SPA : spectrum analyzer, MWVNA: vector network analyzer
        #self._gpib.baud_rate = 57600
        return(self._gpib.query(':INST?'))
    def __del__(self):
        self._gpib.close()
    def idn(self):
        return(self._gpib.query('*IDN?'))
    def avgclear(self):
        self._gpib.write(':SENS:AVER:CLE')
    def avgcount(self, p):
        self._gpib.write(':SENS:AVER:COUN ' + str(p))
    def start(self, p):
        self._gpib.write(':SENS:FREQ:STAR ' + str(p))
        self._gpib.baud_rate = 57600
    def stop(self, p):
        self._gpib.write(':SENS:FREQ:STOP ' + str(p))
        self._gpib.baud_rate = 57600
    def startq(self):
        st = self._gpib.query(':SENS:FREQ:STAR?')
        print(st)
    def sweep(self, sweeptype):
        if sweeptype =='single':
            self._gpib.write(':INITiate:HOLD OFF')
            self._gpib.write(':SENSe:SWEep:TYPE SINGle')
        elif sweeptype == 'cont':
            self._gpib.write(':INITiate:HOLD OFF')
            self._gpib.write(':SENSe:SWEep:TYPE CONTinuous')  # defalut val
    def IFBW(self, p):
        self._gpib.write(':SENS:SWE:IFBW ' + str(p))  # must be one of 100000|50000|20000|10000|5000|2000|1000|500|200|100|50|20|10
        self._gpib.baud_rate = 57600
    def points(self, p):
        self._gpib.write(':SENSe:SWE:POINt ' + str(p))
        self._gpib.baud_rate = 57600
    def minimum(self):
        mini = self._gpib.write(':CALCulate:MARKer 1 :MINimum')
        print(mini)
    def save_s2p(self, filename): #save s2p file
        self._gpib.write(':MMEMory:MSIS INTernal')                  # VERY important to define temp data storage first!
        self._gpib.timeout = 120000                                 # Anritsu is really lazy.. give it sufficient time to respond...:( 
        self._gpib.write(':MMEM:STOR:TRAC 4, \"tempdata\"')
        self._gpib.baud_rate = 57600
        time.sleep(1)
        self._gpib.timeout = 120000
        data = self._gpib.query(':MMEM:DATA? \"tempdata.s2p\"')
        self._gpib.baud_rate = 57600
        self._gpib.timeout = 120000
        time.sleep(2)
        data = data.replace('\r','')
        data = data[data.index('!'):]
        time.sleep(1)
        with open(filename,'w') as f:
            f.write(data[data.index('!'):])
            
    def save_trace(self,filename):
        temp = time.strftime("%d%b%Y%H%M%S", time.localtime())
        self._gpib.write(':MMEMory:MSIS INTernal')
        self._gpib.timeout = 120000     
        self._gpib.write(':MMEM:STOR:TRAC 0,\"'+str(temp)+'\"')
        self._gpib.baud_rate = 57600
        self._gpib.timeout = 120000   
        data = self._gpib.query(':MMEMory:DATA? \"'+str(temp)+'.spa\"')
        self._gpib.baud_rate = 57600
        self._gpib.timeout = 120000
        data = data.replace('\r','')
        data = data.replace('MHz', '')
        data = data.replace('=', ',')
        data = data[data.index('\nP_'):]
        data = data[:data.index('#')]
        with open(filename,'w') as f:
            f.write(data.replace('P_', ''))
            
    def transfer(self, filename, path): # transfer internal memory to pc
        data = self._gpib.query(':MMEM:DATA? \"'+ filename +'\"')
        data = data.replace('\r','')
        data = data[data.index('!'):]
        with open(path+filename,'w') as f:
            f.write(data[data.index('!'):])
            
# =============================================================================
#             
# # In[anritz debug]
# 
# vna = rm.open_resource('TCPIP::'+'192.168.0.105'+'::INSTR')
# 
# # In[data query]
# 
# vna.write(':MMEMory:MSIS INTernal')
# vna.write(':MMEMory:STORe:TRACe 4, \"tempda234234\"')
# vna.baud_rate = 57600
# time.sleep(1)
# 
# data = vna.query(':MMEM:DATA? \"tempda234234.s2p\"')
# vna.baud_rate = 57600
# time.sleep(1)
# 
# # In[dt]
# print(data)
# 
# # In[data query]
# dis = str(time.strftime("%d %b %Y %H:%M:%S", time.localtime()))
# # In[]
# vna.write(':MMEMory:MSIS INTernal')
# vna.write(':MMEMory:STORe:TRACe 0, \"tempda234234 '+dis+'\"')
# vna.baud_rate = 57600
# time.sleep(1)
# 
# # In[]
# vna.timeout = 120000
# data = vna.query(':MMEM:DATA? \"tempda234234 '+dis+'.spa\"')
# vna.baud_rate = 57600
# vna.timeout = 120000
# time.sleep(1)
# 
# # In[dt]
# print(data)
# 
# # In[]
# =============================================================================

        
### Keysight E5071C vector network analyzer
class E5071C:
    def __init__(self, addr):
        self._gpib = rm.open_resource('TCPIP::'+str(addr)+'::INSTR')
    def __del__(self):
        self._gpib.close()
    def idn(self):
        return(self._gpib.query('*IDN?'))
    def save_s2p(self, filename): #save s2p file
        self._gpib.write(':MMEM:STOR \"test.s2p\"')
        self._gpib.timeout = 120000
        data = self._gpib.query(':MMEM:TRAN? \"test.s2p\"') #\"'+filename+'\"')
        self._gpib.timeout = 120000
        data = data[data.index('!'):]
        data = data.replace('\r','')
        with open(filename,'w') as f:
            f.write(data[data.index('!'):])
    def trig(self):
        self._gpib.write(':INIT')
    def avecle(self):
        self._gpib.write('SENS:AVER:CLE')
    def Pow(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SOUR:POW ' + p)
    def cent(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:FREQ:CENT ' + p)
    def marker_max(self):
        self._gpib.write('CALC:MARK:FUNC:TYPE MAX')
        self._gpib.write('CALC:MARK:FUNC:EXEC')
        self._gpib.write('CALC:MARK:SET CENT')
        
    def marker_min(self):
        self._gpib.write('CALC:MARK:FUNC:TYPE MIN')
        self._gpib.write('CALC:MARK:FUNC:EXEC')
        self._gpib.write('CALC:MARK:SET CENT')
        
    def span(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:FREQ:SPAN ' + p)
    def start(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:FREQ:START ' + p)
    def stop(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:FREQ:STOP ' + p)
    def IFBW(self,p):
        #p = np.array2string(p)
        p = numtostr(p)
        self._gpib.write('SENS:BWID ' + p )
    def PowerOff(self):
        self._gpib.write('OUTP 0')
        return
    def PowerOn(self):
        self._gpib.write('OUTP 1')
        return
    def SetPoints(self, x):                                 # Not available in PNA-X (Keisght familiy)
        mystr = '%d' % x
        mystr = 'SENS:SWE:POIN ' + mystr
        self._gpib.write(mystr)
    def maxfreq(self):
        data = self._gpib.query('CALC:MARK:X?')
        return(data)
# =============================================================================
# 
# # In[pnax debug]:
# vna= rm.open_resource('TCPIP::192.168.0.103::INSTR')
# 
# vna.write(':MMEM:STOR \"test.s2p\"')
# vna.timeout = 120000
# data = vna.query(':MMEM:TRAN? \"test.s2p\"')
# vna.timeout = 120000
# 
#     
# # In[]
# =============================================================================
#SR lock-in amplifier

class SR830:
    def __init__(self):
        self._gpib = rm.open_resource('GPIB0::8::INSTR')
    def __del__(self):
        self._gpib.close()
    def mag(self):
        self._gpib.write('OUTP? 3')
        x = self._gpib.read()
        v = np.float_(np.fromstring(x,dtype=float,sep=' '))
        return v
    
### EG&G 5210 Lock-in Amplifier
class EGG5210:
    def __init__(self, addr):
        self._gpib = rm.open_resource('GPIB0::'+str(addr)+'::INSTR')
    def __del__(self):
        self._gpib.close()
    def idn(self):
        return(self._gpib.query("ID"))
    def setvolt(self, v):
        self._gpib.write('OA '+str(v))  # mV
    def setfreq(self, f):
        decimal = int(np.around(np.log10(f),0))
        unit = 4-decimal
        F = int(f*10**(unit))
        BW = str(decimal)
        amp = str(F)
        self._gpib.write('OF '+amp+' '+BW)
    def meas_X(self):
        self._gpib.write('X')
        x_v = self._gpib.read()
        v = np.float_(np.fromstring(x_v ,dtype=float,sep=' '))    # mV
        return v
    def meas_Y(self):
        self._gpib.write('Y')
        y_v = self._gpib.read()
        v = np.float_(np.fromstring(y_v ,dtype=float,sep=' '))  # mV
        return v  
    def meas_XY(self):
        self._gpib.write('XY')
        xy = self._gpib.read()
        v = np.float_(np.fromstring(xy ,dtype=float,sep=' '))  # mV
        return v
    def meas_mag(self):
        self._gpib.write('MAG')
        mag = self._gpib.read()
        v = np.float_(np.fromstring(mag ,dtype=float,sep=' '))  # mV
        return v 
    def meas_phase(self):
        self._gpib.write('PHA')
        phase = self._gpib.read()
        ph = np.float_(np.fromstring(phase ,dtype=float,sep=' '))  # millidegree
        return ph/1e3
    
        
    
### xilinx RFSoC
'''
nameserver (type on PUTTY): 
    PYRO_SERIALIZERS_ACCEPTED=pickle PYRO_PICKLE_VERSION=4 pyro4-ns -n 192.168.0.102 -p 8888

nameserver list (type on local client): 
    pyro4-nsc -v list
    
server (type on PUTTY): 
    python xilinx/jupyter_notbooks/qick/pyro4/server2.py 
    or open http://192.168.0.102:9090/notebooks/qick/pyro4/Remote_operation_with_Pyro4.ipynb 
    and run server code
 
'''

class RFSoC():
    def __init__(self):
        Pyro4.config.SERIALIZER = "pickle"
        Pyro4.config.PICKLE_PROTOCOL_VERSION=4
        
        self.ns_host = "192.168.0.102"
        self.ns_port = 8888
        self.server_name = "myqick"
        
        self.ns = Pyro4.locateNS(host=self.ns_host, port=self.ns_port)


        for k,v in self.ns.list().items():
            print(k,v)

        
        self.soc = Pyro4.Proxy(self.ns.lookup(self.server_name))
        self.soccfg = QickConfig(self.soc.get_cfg())
        print(self.soccfg)
        
    def socconfig(self):
        cfg = Pyro4.Proxy(self.ns.lookup(self.server_name))
        soc = QickConfig(self.soc.get_cfg())
        return cfg, soc

from qick import *

class SinglePulse(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   
        res_ch = cfg["res_ch"]

        # set the nyquist zone
        self.declare_gen(ch=cfg["res_ch"], nqz=1)
        
        # configure the readout lengths and downconversion frequencies (ensuring it is an available DAC frequency)
        self.declare_readout(ch=cfg["ro_ch"], length=self.cfg["readout_length"],
                                 freq=self.cfg["readout_freq"], gen_ch=cfg["res_ch"])

        # convert frequency to DAC frequency (ensuring it is an available ADC frequency)
        freq = self.freq2reg(cfg["pulse_freq"],gen_ch=res_ch, ro_ch=cfg["ro_ch"])
        phase = self.deg2reg(cfg["res_phase"], gen_ch=res_ch)
        gain = cfg["pulse_gain"]
        self.default_pulse_registers(ch=res_ch, freq=freq, phase=phase, gain=gain)

        style=self.cfg["pulse_style"]

        if style in ["flat_top","arb"]:
            sigma = cfg["sigma"]
            self.add_gauss(ch=res_ch, name="measure", sigma=sigma, length=sigma*5)
            
        if style == "const":
            self.set_pulse_registers(ch=res_ch, style=style, length=cfg["length"], mode=cfg["mode"])
        elif style == "flat_top":
            # The first half of the waveform ramps up the pulse, the second half ramps down the pulse
            self.set_pulse_registers(ch=res_ch, style=style, waveform="measure", length=cfg["length"])
        elif style == "arb":
            self.set_pulse_registers(ch=res_ch, style=style, waveform="measure")
        
        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        # fire the pulse
        # trigger all declared ADCs
        # pulse PMOD0_0 for a scope trigger
        # pause the tProc until readout is done
        # increment the time counter to give some time before the next measurement
        # (the syncdelay also lets the tProc get back ahead of the clock)
        
        self.measure(pulse_ch=self.cfg["res_ch"], 
                     adcs=self.ro_chs,
                     pins=[0], 
                     adc_trig_offset=self.cfg["adc_trig_offset"],
                     wait=True,
                     syncdelay=self.us2cycles(self.cfg["relax_delay"]))
        

class LoopbackProgram(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   
#         self.declare_gen(ch=cfg["jpa_ch"], nqz=1) #JPA
        self.declare_gen(ch=cfg["res_ch"], nqz=1) #Readout
#         self.declare_gen(ch=cfg["qubit_ch"], nqz=2) #Qubit
#         self.declare_gen(ch=cfg["storage_ch"], nqz=2) #Storage

        for ch in [0,1]: #configure the readout lengths and downconversion frequencies
            self.declare_readout(ch=ch, 
                                 length=cfg["readout_length"],
                                 freq=cfg["frequency"], 
                                 gen_ch=cfg["res_ch"])

        freq=self.freq2reg(cfg["frequency"], 
                           gen_ch=cfg["res_ch"], ro_ch=0)  # convert frequency to dac frequency (ensuring it is an available adc frequency)
        self.set_pulse_registers(ch=cfg["res_ch"], 
                                 style="const", 
                                 freq=freq, 
                                 phase=0, 
                                 gain=cfg["pulse_gain"],
                                length=cfg["pulse_length"])
        
        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        cfg=self.cfg   
        self.measure(pulse_ch=cfg["res_ch"], 
             adcs=[0,1],
             adc_trig_offset=cfg["adc_trig_offset"],
             t=0,
             wait=True,
             syncdelay=self.us2cycles(cfg["relax_delay"]))        



class SingleToneSpectroscopyProgram(AveragerProgram):
    def initialize(self):
        
        cfg=self.cfg   
        self.declare_gen(ch=cfg["res_ch"], nqz=1) #Readout
        
        #for ch in [0,1]: #configure the readout lengths and downconversion frequencies
        
        self.declare_readout(ch=cfg["ro_ch"], 
                             length=cfg["readout_length"],
                                 freq=cfg["frequency"], 
                                 gen_ch=cfg["res_ch"])
        
        freq=self.freq2reg(cfg["frequency"], 
                           gen_ch=cfg["res_ch"], 
                           ro_ch=cfg["ro_ch"])  # convert frequency to dac frequency (ensuring it is an available adc frequency)
        
        self.set_pulse_registers(ch=cfg["res_ch"], 
                                 style="const", 
                                 freq=freq, 
                                 phase=cfg["res_phase"], 
                                 gain=cfg["res_gain_start"],
                                length=cfg["readout_length"])

        self.synci(200)  # give processor some time to configure pulses

    def body(self):
        self.measure(pulse_ch=self.cfg["res_ch"], 
             #adcs=[0,1],
             adcs=self.ro_chs,
             adc_trig_offset=self.cfg["adc_trig_offset"],
             wait=True,     
             syncdelay=self.us2cycles(self.cfg["relax_delay"]))  
    

    

        

            

