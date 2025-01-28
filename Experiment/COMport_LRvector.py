# In[0] : necessary packages
import numpy as np

import serial

import json
import time
import sys


# In[test @ heliox pc]

com = serial.Serial('COM49', baudrate=115200, timeout=None)

com.flushInput()
com.flushOutput()

# In[]

#com.write(str.encode('?'))
com.write(bytes('?\n', 'utf-8'))

# In[]

data = com.readline().rstrip()
print(data)

com.flushInput()
com.flushOutput()


# In[]

com.write(bytes('RF 1\n', 'utf-8'))

# =============================================================================
# com.write(bytes('DD 1\n', 'utf-8'))
# 
# com.write(bytes("WF 1 0 0\n", 'utf-8'))
# =============================================================================

# In[]

com.write(bytes('RW 1\n', 'utf-8'))

# In[]:
    
'''
"IC dac Vret Vmin Vmax Vthp Vthn t1 t2 t3 tramp repeats AIN" 
'''

com.write(bytes('IC 1 0 -1 1 12 -12 0.01 0.01 0.001 1.0 1 1\n', 'utf-8'))

# In[]

com.write(bytes('RC 1\n', 'utf-8'))
data
# In[]

com.flushInput()
com.flushOutput()


com.write(bytes('IS 1 2\n', 'utf-8'))

# In[]:
    
while True:
    com.rts = False
    com.dtr = False
    bytesToRead = com.inWaiting()
    data=com.read(bytesToRead)
    print(data)

# In[]

com.close()

# In[DAC]

'''
DAC test
'''

com = serial.Serial('COM4', baudrate=115200, timeout=None)

com.flushInput()
com.flushOutput()

# In[]

com.write(str.encode('?\n'))

# In[]

data = com.readline().rstrip()
#print(1)
print(data)

# In[]

com.write(str.encode('RF 1\n'))
time.sleep(0.05)
data = com.readline().rstrip()
print(data)

# In[]

data = com.readline().rstrip()
#print(1)
print(data)
    
# In[]

data=com.read(size=1)

print(data)

# In[]

com.flushInput()
com.flushOutput()

com.write(str.encode('DD 1\n'))

# In[]



com.close()

# In[magnet]

'''
magnet test
'''

com = serial.Serial('COM8', baudrate=115200, timeout=0.1)

com.flushInput()
com.flushOutput()

# In[]
seg =1
com.write(str.encode(f'STATE? \n'))

# In[]

data = com.readline().rstrip()
#print(1)
print(data)

# In[]

com.write(str.encode('RF 1\n'))
time.sleep(0.05)
data = com.readline().rstrip()
print(data)

# In[]

data = com.readline().rstrip()
#print(1)
print(data)
    
# In[]

data=com.read(size=1)

print(data)

# In[]

com.flushInput()
com.flushOutput()

com.write(str.encode('DD 1\n'))

# In[]



com.close()



