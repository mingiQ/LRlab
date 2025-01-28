# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import SerialInstrument



## MicroLambda bandpass filter
class MLBPF(SerialInstrument):
    '''
        /*    ------ communication commands ---------
     *     serial @ 115200, line end character: \n
     *    "?" or "*IDN?"             ->  devstamp         - equipment ID  
     *    "FF val"                   -> "FF"              - set single dac value (float)
     *    "SF fmin fmax"             -> "SF"              - save fmin and fmax values to EEPROM
     *    "RF"                       -> "RF"              - read fmin and fmax values from EEPROM
     *    
     *    "TS dac dval"              -> "TS"              - comm test: write xAAAAA to dac dval times
     *    "TX pin dval"              -> "TX"              - pin test: write HIGH/LOW to pin dval times
     *    "TT pin"                   -> "TT"              - toggle pin
    */
    '''
    def __init__(self, name="ML_BPF", address='COM10', enabled=True, timeout=0.1, baudrate=115200):
        SerialInstrument.__init__(self, name, address, enabled, timeout, baudrate)
        self.flush_input()
        self.flush_output()
        
    def idn(self):
        return self.query('?')
    
    def freq(self, f):
        return self.query(f'FF {f}')
    
    def debug(self, f):
        return self.query(f'FH {f}')
    
# =============================================================================
# 
# # In[]
# 
# if __name__ == '__main__':
#     bpf = MLBPF('BPF', address='COM3')
#     print(bpf.idn())
#     print('taking data')
#     
# # In[test]
# 
# bp = MLBPF(address='COM3')
# 
# print(bp.idn())
# 
# 
# =============================================================================
