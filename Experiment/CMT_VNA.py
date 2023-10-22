# In[0]:

import pyvisa
import os
import numpy as np

# In[1]

rm = pyvisa.ResourceManager('C:/Windows/System32/visa32.dll')

# In[connet cmt vna]

cmt = rm.open_resource('TCPIP0::127.0.0.1::5025::SOCKET')

# In[]

cmt.write(':SENS1:FREQ:STAR 3 GHz \n')
cmt.write(':SENS1:FREQ:STOP 5.5 GHz \n')

# In[]

cmt.write(':SENS1:SWE:POIN 20000 \n')
cmt.write(':SOUR:POW 0 dBm \n')

# In[]:
    
cmt.write(':OUTP 1')

# In[]
cmt.write(':INIT  \n')

# In[]
cmt.write(':TRIG:SOUR BUS \n')

cmt.write(':TRIG:SEQ:SING \n')


# In[save]

cmt.write(':MMEM:STOR:SNP \"test2.s2p\" \n')
cmt.timeout = 120000

# In[trasfer tenp file to experiment folder]
expt_dir = 'Z:/Mingi Kim/Resonators/230813_Leiden_Mount/CM/'

data = open("C:/Program Files/S2VNA/FixtureSim/test2.s2p", 'r')

with open("C:/Program Files/S2VNA/FixtureSim/test2.s2p", 'r') as source, open(expt_dir+'test.s2p', 'w') as target:
    for line in source:
        target.write(line)



# In[transfer data -- not working?]
dat = cmt.query(':MMEM:TRAN? \"C:\Program Files\S2VNA\FixtureSim\test2.s2p\" \n')
cmt.timeout = 120000

# In[]

