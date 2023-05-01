# In[0] : necessary pkgs

import socket as socket
import numpy as np
import pandas as pd
import os
import time

# In[0.5] : leiden log

Llog = 'C:/CF81-logs/'


# In[1] : 
    
cmdlist = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

# server code

HOST = '10.164.26.139'

PORT = 8888

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
print('socket created')

time.sleep(1)
print('wait until socket activated ... ')

for sec in range(5):
    print(5-sec)
    time.sleep(1)
    
try:
    print("Binding socket to port: " + str(PORT))
    s.bind((HOST, PORT))
    
except socket.error as msg:
    print('Bind failed:' +str(msg))
    print("more time is required ... ")
    wait = input('Wait ten more seconds? (Y/N)')
    
    if str(wait) == 'Y':
        time.sleep(10)
        s.bind((HOST, PORT))
        
    elif str(wait) == 'N':
        print('BIND ERR: ' + str(msg))
        
s.listen(5)

print('Socket awaiting messages')
conn, addr = s.accept()
print('Connected with ' + str(addr))
msg = 'CONNECTED!'
conn.send(msg.encode())

# sever while loop

while True:
    
    cmd = conn.recv(1024)
    
    print('I sent a msg back in response to: ' + cmd.decode())
    reply = ''
    
    "read temperature data"
    
    if cmd.decode() in cmdlist:
        files = os.listdir(Llog)[-10::]  # most recent 10 TC logs
        ind = int(cmd.decode())
        tempdata = np.loadtxt(fname=Llog+files[ind], skiprows=4, usecols=[ind+13])[-1]
        reply = str(tempdata)

    elif cmd.decode() == 'quit':
        conn.close()

    elif cmd.decode() == 'exit':
        conn.close()
        break

    else:
        reply = 'Unknown command'

    conn.send(reply.encode())

conn.close()
    
# In[debugs]

files = os.listdir(Llog)[-10::]
temp = np.loadtxt(fname=Llog+files[9], skiprows=4, usecols=[22])
print(temp)
    

