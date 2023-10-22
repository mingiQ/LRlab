# In[0] : necessary packages
import numpy as np
import pyvisa
import serial
import socket
import json
import time
import sys


# In[]

com = serial.Serial('COM11', baudrate=115200, timeout=None)

com.flushInput()
com.flushOutput()

# In[]

#com.write(str.encode('?'))
com.write(bytes('?', 'utf-8'))

# In[]

data = com.readline().rstrip()
print(data)

com.flushInput()
com.flushOutput()

# In[]:
    
while True:
    com.rts = False
    com.dtr = False
    bytesToRead = com.inWaiting()
    data=com.read(bytesToRead)
    print(data)

# In[]

com.close()

# In[]

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

