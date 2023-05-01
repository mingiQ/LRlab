# In[1]:
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from scipy.optimize import curve_fit
from numpy import log10, absolute, sqrt

# In[2022/11/30 : directory]
path = 'C:/Users/20141/OneDrive - purdue.edu/PhD_research/resonators/1129_Leiden_mount/1130_biased resonator/'
data_folder = path + 'S21data/'

# In[data]

def data_extract(Path, file):
    dat = np.loadtxt(Path+file, skiprows=23)
    freq = dat[:,0]
    s21 = dat[:,3]
    phase = dat[:,4]
    return (freq, s21, phase)

# In[plot]
filelist = os.listdir(data_folder)
datlist_50dBaten = filelist[:11]
datlist_30dBaten_20kO = filelist[12:22]
datlist_30dBaten_100kO = filelist[23:33]

for files in filelist:
    data = data_extract(data_folder, files)
    plt.subplot(211)
    plt.plot(data[0], data[1])
    plt.xlabel('Frequency(GHz)')
    plt.ylabel('S21(dB)')
    
    plt.subplot(212, projection='polar')
    plt.plot(data[2]*np.pi/180, 10**(data[1]/20) )
    plt.grid(True)

    
    plt.tight_layout()
    


# In[correlation fucntion between spectra]

dat1 = data_extract(data_folder, filelist[22])  # reference signal
dat2 = data_extract(data_folder, filelist[23])  # test signal
freq_renorm = np.linspace(dat1[0][0], dat1[0][-1], 4001)

slices = len(dat1[0])

'''
say, integration region to be 4.4325GHz to 4.4331GHz (linewidth ~ 200kHz) larger than cavity linewidth

I matched frequency step with 262.5Hz / total dataset starts from 4.4322GHz to 4.43325GHz, 4001 slices, (index runs from 0 to 4000)
'''
step = 2.625e-7


'''
cross-correlation function at f=-1000*step

'''
fAB=[]
SAB=[]

fAB.append(-1000*step)

'''
SA --> summation ranges over 1000 ~ 3000
SB --> summation ranges over 2000 ~ 4000

'''
plt.plot(freq_renorm[1000:3000], (10**(dat1[1]/20))[1000:3000])  # ref
plt.plot((freq_renorm-1000*step)[2000:4000], (10**(dat1[1]/20))[2000:4000]) # input

s=0
for i in range(len(freq_renorm[1000:3000])):
   s+= (10**(dat1[1]/20))[1000:3000][i]*(10**(dat2[1]/20))[2000:4000][i]*step
   
SAB.append(s)

# In[cross-correlation function]
scan_range = 1000


sum_range = np.arange(-scan_range,scan_range+1,1)
step = 2.625e-7  # 0.26kHz
fAB=[]
SAB=[]

for i in range(len(sum_range)):
    order = sum_range[i]
    fAB.append(order*step)
    s = 0
    finit = 2000-scan_range
    fend = 2000+scan_range
    for ii in range(len(freq_renorm[finit:fend])):
        s+= (10**(dat1[1]/20))[finit:fend][ii]*(10**(dat1[1]/20))[finit-order:fend-order][ii]*step  # let's check self-correlation
        
    SAB.append(s)

# In[cross corrleation plot]
plt.plot(np.array(fAB)*1e6, np.log(np.array(SAB)), '.')
plt.xlabel('frequency shift from zero-bias cavity(kHz)')
plt.ylabel(r'Cross-correlation $S_{II_0}(\Delta f)$')

# In[cross-correlation function]

correlat_path = 'C:/Users/20141/OneDrive - purdue.edu/PhD_research/resonators/1129_Leiden_mount/1130_biased resonator/correlaton function/'

def cross_correlate(dat1, dat2, scan_range, filename):

    sum_range = np.arange(-scan_range,scan_range+1,1)
    freq_renorm = np.linspace(dat1[0][0], dat1[0][-1], 4001)
    step = 2.625e-7  # 0.26kHz
    fAB=[]
    SAB=[]
    
    for i in range(len(sum_range)):
        order = sum_range[i]
        fAB.append(order*step)
        s = 0
        finit = 2000-scan_range
        fend = 2000+scan_range
        for ii in range(len(freq_renorm[finit:fend])):
            s+= (10**(dat1[1]/20))[finit:fend][ii]*(10**(dat2[1]/20))[finit-order:fend-order][ii]*step  # let's check self-correlation
            
        SAB.append(s)
        
    data = np.vstack((fAB, SAB))
    np.savetxt(correlat_path+filename+'.csv', data)
    
    plt.plot(np.array(fAB)*1e6, np.log(np.array(SAB)), '.')
    plt.xlabel('frequency shift from zero-bias cavity(kHz)')
    plt.ylabel(r'Cross-correlation $S_{II_0}(\Delta f)$')

def square_correlate(dat1, dat2, scan_range, filename):

    sum_range = np.arange(-scan_range,scan_range+1,1)
    freq_renorm = np.linspace(dat1[0][0], dat1[0][-1], 4001)
    step = 2.625e-7  # 0.26kHz
    fAB=[]
    SAB=[]
    
    for i in range(len(sum_range)):
        order = sum_range[i]
        fAB.append(order*step)
        s = 0
        finit = 2000-scan_range
        fend = 2000+scan_range
        for ii in range(len(freq_renorm[finit:fend])):
            s+= ((10**(dat1[1]/20))[finit:fend][ii]-(10**(dat2[1]/20))[finit-order:fend-order][ii])**2*step  # let's check self-correlation
            
        SAB.append(s)
        
    data = np.vstack((fAB, SAB))
    np.savetxt(correlat_path+filename+'.csv', data)
    
    plt.plot(np.array(fAB)*1e6, np.log(np.array(SAB)), '.')
    plt.xlabel('frequency shift from zero-bias cavity(kHz)')
    plt.ylabel(r'Square-correlation $log\,S_{II_0}(\Delta f)$')
    

# In[cross corrleation plot]
cross_correlate(data_extract(data_folder, filelist[22]), data_extract(data_folder, filelist[23]), 1000, 'I0vsI500uA')

# In[cross correlation plot]
for i in range(9):
    cross_correlate(data_extract(data_folder, filelist[22]), data_extract(data_folder, filelist[24+i]), 1000, 'I0vsI'+str(50*(i+1))+'uA')
    print('finish {}th data'.format(i))
    
# In[squqre correlation plot]

for i in range(11):
    square_correlate(data_extract(data_folder, filelist[22]), data_extract(data_folder, filelist[22+i]), 1000, 'square_cor_I0vsI'+str(50*(i))+'uA')
    print('finish {}th data'.format(i))    
    
# In[saved data]
datalist = os.listdir(correlat_path)

plt.figure(figsize=(10,5))
for data in datalist:
    dat = np.loadtxt(correlat_path+data)
    plt.plot(dat[0]*1e6, dat[1], '.', label=data)
    #plt.plot(np.array(fAB)*1e6, np.log(np.array(SAB)), '.')
    plt.xlabel('frequency shift from zero-bias cavity(kHz)')
    plt.ylabel(r'Cross-correlation $S_{II_0}(\Delta f)$')
    plt.legend(loc='upper center')
    
    
# In[density plot]
vlist = np.arange(1,11,1)
a = 400

frequency = (dat[0][1000-a:1001+a])*1e6

X,Y = np.meshgrid(frequency, vlist)

Z = np.zeros((len(vlist), len(frequency)))

for v in range(len(vlist)):
    file='I0vsI'+str(50*(v+1))+'uA.csv'
    dat = np.loadtxt(correlat_path+file)
    Z[v,:] = dat[1][1000-a:1001+a]
    
# In[5] Plot density plot
plt.pcolormesh(X,Y/20,np.log(Z))
plt.xlabel('Frequency(kHz)')
plt.ylabel('bias current(mA)')

cbar = plt.colorbar()
cbar.set_label(r'Cross-correlation $log\,S_{II_0}(\Delta f)$', rotation=270)
plt.tight_layout()

# In[density plot]
vlist = np.arange(1,11,1)
a = 80
dat0 = np.loadtxt(correlat_path+'square_cor_I0vsI0uA.csv')
frequency = (dat0[0][1000-a:1001+a])*1e6

X,Y = np.meshgrid(frequency, vlist)

Z = np.zeros((len(vlist), len(frequency)))

for v in range(len(vlist)):
    file='square_cor_I0vsI'+str(50*(v+1))+'uA.csv'
    dat = np.loadtxt(correlat_path+file)
    Z[v,:] = dat[1][1000-a:1001+a]
    
# In[5] Plot density plot
plt.pcolormesh(X,Y/20,np.log(Z))
plt.xlabel('Frequency(kHz)')
plt.ylabel('bias current(mA)')

cbar = plt.colorbar()
cbar.set_label(r'Cross-correlation $log\,S_{II_0}(\Delta f)$', rotation=270)
plt.tight_layout()



