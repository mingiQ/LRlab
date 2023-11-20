# In[0]:
    
__author__ = 'Mingi'

import sys
import time
import os
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments import InstrumentManager, mlBPF, TestInstrument



# In[1]

im = InstrumentManager()



bpf = mlBPF.MLBPF("BPF", address='COM3', timeout=0.1)
test = TestInstrument.EchoInstrument('test', address='addr')

# In[]

bpf.idn()