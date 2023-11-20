# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import VisaInstrument


## Anritsu MS2038C VNA master  ''' TCPIP[board]::host address[::LAN device name][::INSTR] '''


class MS2038(VisaInstrument):
    
    def __init__(self, name='MS2038c', address='TCPIP::192.168.0.105::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=1
        self.recv_length=65536
        self.term_char=''
   
    def idn(self):
         return self.instrument.query('*IDN?')
         
    def mode(self, mod): 
        self.instrument.write(':INST[:SEL] ' + str(mod))  # SPA : spectrum analyzer, MWVNA: vector network analyzer
        #self._gpib.baud_rate = 57600
        return(self.instrument.query(':INST?'))
    
    def start(self, p):
        self.instrument.write(':SENS:FREQ:STAR ' + str(p))
        self.instrument.baud_rate = 57600
    def stop(self, p):
        self.instrument.write(':SENS:FREQ:STOP ' + str(p))
        self.instrument.baud_rate = 57600
    def startq(self):
        st = self.instrument.query(':SENS:FREQ:STAR?')
        print(st)
    def IFBW(self, p):
        self.instrument.write(':SENS:SWE:IFBW ' + str(p))  # must be one of 100000|50000|20000|10000|5000|2000|1000|500|200|100|50|20|10
        self.instrument.baud_rate = 57600
    def points(self, p):
        self.instrument.write(':SENSe:SWE:POINt ' + str(p))
        self.instrument.baud_rate = 57600
    
    def avgclear(self):
        self.instrument.write(':SENS:AVER:CLE')
    def avgcount(self, p):
        self.instrument.write(':SENS:AVER:COUN ' + str(p))
   
    def sweep(self, sweeptype):
        if sweeptype =='single':
            self.instrument.write(':INITiate:HOLD OFF')
            self.instrument.write(':SENSe:SWEep:TYPE SINGle')
        elif sweeptype == 'cont':
            self.instrument.write(':INITiate:HOLD OFF')
            self.instrument.write(':SENSe:SWEep:TYPE CONTinuous')  # defalut val
    
    def minimum(self):
        mini = self.instrument.write(':CALCulate:MARKer 1 :MINimum')
        print(mini)
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
    
# In[]

vna = MS2038()   
print(vna.idn())    
   
