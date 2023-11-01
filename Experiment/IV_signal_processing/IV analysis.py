# In[0]: 
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import pandas as pd

sys.path.append('//lrgrouphd2/usr/general/LRlabcode/LRlab/Experiment')
#sys.path.append('Z:/general/LRlabcode/LRlab/Experiment/Instruments')
#sys.path.append('//lrgrouphd2/usr/Mingi Kim/python_codes/Experiment')

from data_analysis_codes import *
from fit_codes import *
from scipy.constants import h, e

# In[1]:
    
wdir = r'\\LRGROUPHD2\usr\Mingi Kim\0_P2_ABS_spectrometer\Shapiro step\I-V test for dV-dI\\'

ivdata = os.listdir(wdir)

iv1000 = [x for x in ivdata if '1000pts' in x]
iv5000 = [x for x in ivdata if '5000pts' in x]
iv10000 = [x for x in ivdata if '10000pts' in x]
iv_fil = [x for x in ivdata if 'filtered' in x]

# In[import data]

for i in range(len(ivdata)):
    if 'filtered' in ivdata[i]:
        sample = np.loadtxt(wdir+ivdata[i], skiprows=2, delimiter=',')
    else:
        sample = np.loadtxt(wdir+ivdata[i], delimiter=',')
        
    plt.figure(figsize=(5,5))
    plt.plot(sample[:,0], sample[:,1], label=ivdata[i])
    plt.legend()
    
    
# In[]
ivlist = iv10000
for i in range(len(ivlist)):
    sample = np.loadtxt(wdir+ivlist[i], delimiter=',')
    print(i)
        
    plt.figure(figsize=(5,5))
    plt.plot(sample[:,0], sample[:,1], label=ivlist[i])
    plt.plot(sample[:,0], sample[:,2], label=ivlist[i])
    plt.legend()
    
# In[]
from scipy.fft import fft, fftfreq, fftshift
from scipy.signal import lfilter, decimate, filtfilt

iv = iv10000[-1]

sample = np.loadtxt(wdir+iv, delimiter=',')
I = sample[:,0]
V = sample[:,1]
Vf = fft(V)
If = fftfreq(len(I), np.abs(I[1]-I[2]))

Vfs = fftshift(Vf)
Ifs = fftshift(If)
plt.plot(Ifs, Vfs,'.')
plt.yscale('log')


# In[apply infinite impulse resonsoe filter]


#apply the filter and plot the denoised data
n = 15  # the larger n is, the smoother curve will be
b = [1.0 / n] * n
a = 1
yy = lfilter(b,a,V)
plt.plot(I,V,alpha=0.4)
plt.plot(I, yy, 'b.', markersize=0.01)

# In[decimation]
downsampling_factor = 50

Vdem = decimate(V, downsampling_factor)
Inew = np.linspace(I[0], I[-1], len(Vdem))

plt.plot(Inew, Vdem, '.')


# In[]

dvdi = (Vdem[1:]-Vdem[:-1])/(Inew[1:]-Inew[:-1])
plt.plot(dvdi,'.-')

# In[tesing IIR --> decimation]

iv = iv5000[3]

sample = np.loadtxt(wdir+iv, delimiter=',')
I_raw = sample[:,0]
V_raw = sample[:,1]
fs = 30

b, a = scipy.signal.iirfilter(4, Wn=2.5, fs=fs, btype="low", ftype="butter")
print(b, a, sep="\n")

V_iir = filtfilt(b,a,V_raw)

downsampling_factor = int(len(V_iir)/50)

Vdem = decimate(V_iir, downsampling_factor)
Inew = np.linspace(I[0], I[-1], len(Vdem))

dvdi = (Vdem[1:]-Vdem[:-1])/(Inew[1:]-Inew[:-1])
ilist = 0.5*(Inew[1:]+Inew[:-1])

plt.figure(figsize=(5,10))
plt.subplot(211)
plt.plot(Inew, Vdem, '.')
plt.plot(I_raw, V_iir, 'g-', alpha=0.4)
plt.plot(I_raw, V_raw, 'r-', alpha=0.4)

plt.subplot(212)
plt.plot(ilist, dvdi, '.-')

# In[ decimation]

for i in range(6)
    iv = iv5000[i]
    
    sample = np.loadtxt(wdir+iv, delimiter=',')
    I_raw = sample[:,0]
    V_raw = sample[:,1]
    
    V_iir = V_raw
    downsampling_factor = int(len(V_iir)/50)
    
    Vdem = decimate(x=V_iir, q=downsampling_factor, n=3, ftype='iir' ,zero_phase=True)
    Inew = np.linspace(I[0], I[-1], len(Vdem))
    
    dvdi = (Vdem[1:]-Vdem[:-1])/(Inew[1:]-Inew[:-1])
    ilist = 0.5*(Inew[1:]+Inew[:-1])
    
    plt.figure(figsize=(5,10))
    plt.subplot(211)
    plt.plot(Inew, Vdem, '.', label=iv)
    plt.plot(I_raw, V_iir, 'g-', alpha=0.4)
    plt.plot(I_raw, V_raw, 'r-', alpha=0.4)
    plt.legend()
    
    plt.subplot(212)
    plt.plot(ilist, dvdi, '.-')


# In[]
import seaborn as sns
sns.set_theme(style="whitegrid")

# Load the example exercise dataset
exercise = sns.load_dataset("exercise")

# Draw a pointplot to show pulse as a function of three categorical factors
g = sns.catplot(
    data=exercise, x="time", y="pulse", hue="kind", col="diet",
    capsize=.2, palette="YlGnBu_d", errorbar="se",
    kind="point", height=6, aspect=.75,
)
g.despine(left=True)
