# In[0]:
import numpy as np
import pandas as pd
import sys
import os
import time
import matplotlib.pyplot as plt

sys.path.append('Z:/Mingi Kim/python_codes')
import LAB_v0


# In[import temperature controller]

T = LAB_v0.Thermo_leiden()

# In[read 3k pot (readT ch.0)]

now = time.time()
tlist = []
Tlist = []

tlist.append(0)
Tlist.append(T.readT(0))
print(Tlist[-1])

for i in range(10):
    time.sleep(30)
    interval = time.time()-now
    tlist.append(interval)
    Tlist.append(T.readT(0))
    print(Tlist[-1])
    
# In[plot]
plt.plot(tlist, Tlist, '.-')
plt.xlabel('Time(s)')
plt.ylabel('temperature(mK)')

# In[read cold plate (readT ch.2)]

now = time.time()
tlist = []
Tlist = []

tlist.append(0)
Tlist.append(T.readT(2))
print(Tlist[-1])

for i in range(10):
    time.sleep(30)
    interval = time.time()-now
    tlist.append(interval)
    Tlist.append(T.readT(2))
    print(Tlist[-1])
    
# In[plot]
plt.plot(tlist, Tlist, '.-')
plt.xlabel('Time(s)')
plt.ylabel('temperature(mK)')

# In[]

