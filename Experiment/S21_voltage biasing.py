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

# In[data folder]
os.makedirs(data_folder)

# In[measurement function]

def S21trans(start, stop, IF, avg, meas_t, path, name):
    vna.avgcount(avg)                       # average count 
    vna.IFBW(IF)                         # IFBW 2kHz
    vna.start(start)
    vna.stop(stop)
    
    filename = path+name
    
    vna.avgclear()
    time.sleep(meas_t)
    vna.save_s2p(filename)
# In[Anritz]

'''
VNA setting 
'''

vna = LAB_v0.MS2038('192.168.0.105') 

# In[meas]
'''
0V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias0V.s2p')

# In[meas]
'''
1V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias1V.s2p')

# In[meas]
'''
2V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias2V.s2p')

# In[meas]
'''
3V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias3V.s2p')

# In[meas]
'''
4V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias4V.s2p')
# In[meas]
'''
5V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias5V.s2p')

# In[meas]
'''
6V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias6V.s2p')

# In[meas]
'''
7V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias7V.s2p')

# In[meas]
'''
8V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias8V.s2p')

# In[meas]
'''
9V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias9V.s2p')

# In[meas]
'''
10V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 60, data_folder, '10kO_vias10V.s2p')



# In[data]

def data_extract(Path, file):
    dat = np.loadtxt(Path+file, skiprows=23)
    freq = dat[:,0]
    s21 = dat[:,3]
    phase = dat[:,4]
    return (freq, s21, phase)

# In[plot]
filelist = os.listdir(data_folder)
datlist_50dBaten = []
datlist_30dBaten_20kO = []
datlist_30dBaten_100kO = []

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

# In[plot]
filelist = os.listdir(data_folder)
datlist = []
plt.figure(figsize=(5,10))
for files in filelist:
    data = data_extract(data_folder, files)
    plt.subplot(211, projection='polar')
    plt.plot(data[2]*np.pi/180, 10**(data[1]/20) )
    plt.grid(True)
    
    plt.subplot(212)
    plt.plot(data[0], data[2])
    plt.xlabel('Frequency(GHz)')
    plt.ylabel('Phase(DEG)')
    


# In[correlation fucntion between spectra]

dat1 = data_extract(data_folder, filelist[22])  # reference signal
dat2 = data_extract(data_folder, filelist[23])  # test signal
freq_renorm = np.linspace(dat1[0][0], dat1[0][-1], 4001)

#plt.plot(dat1[0], 10**(dat1[1]/20))
#plt.plot(dat1[0], 10**(dat2[1]/20))

#plt.plot(freq_renorm, 10**(dat1[1]/20))
#plt.plot(freq_renorm, 10**(dat2[1]/20))

slices = len(dat1[0])

'''
say, integration region to be 4.4325GHz to 4.4331GHz (linewidth ~ 200kHz) larger than cavity linewidth

I matched frequency step with 262.5Hz / total dataset starts from 4.4322GHz to 4.43325GHz, 4001 slices, (index runs from 0 to 4000)
'''
step = 2.625e-7
# =============================================================================
# 
# df= []
# 
# for i in range(slices-1):
#     df.append(freq_renorm[i+1]-freq_renorm[i])
# 
# =============================================================================

'''
cross-correlation function at f=-1000*step

'''
fAB=[]
SAB=[]

fAB.append(-1000*step)
#plt.plot(freq_renorm, 10**(dat1[1]/20))  # ref
#plt.plot(freq_renorm+1000*step, 10**(dat2[1]/20))
#plt.plot(freq_renorm-1000*step, 10**(dat2[1]/20)) # input

# SB[f+DF] = freq_renorm-1000*step, 10**(dat2[1]/20)
'''
SA --> summation ranges over 1000 ~ 3000
SB --> summation ranges over 2000 ~ 4000

'''
plt.plot(freq_renorm[1000:3000], (10**(dat1[1]/20))[1000:3000])  # ref
#plt.plot(freq_renorm+1000*step, 10**(dat2[1]/20))
plt.plot((freq_renorm-1000*step)[2000:4000], (10**(dat2[1]/20))[2000:4000]) # input

s=0
for i in range(len(freq_renorm[1000:3000])):
   s+= (10**(dat1[1]/20))[1000:3000][i]*(10**(dat2[1]/20))[2000:4000][i]*step
   
SAB.append(s)

# In[cross-correlation function]
sum_range = np.arange(-1000,1001,1)
step = 2.625e-7
fAB=[]
SAB=[]

for i in range(len(sum_range)):
    order = sum_range[i]
    fAB.append(order*step)
    s = 0
    for ii in range(len(freq_renorm[1000:3000])):
        s+= (10**(dat1[1]/20))[1000:3000][ii]*(10**(dat2[1]/20))[1000-order:3000-order][ii]*step
        
    SAB.append(s)

# In[cross corrleation plot]
correlat_path = 'C:/Users/20141/OneDrive - purdue.edu/PhD_research/resonators/1129_Leiden_mount/1130_biased resonator/correlaton function/'

plt.plot(np.array(fAB)*1e6, SAB, '.')
plt.xlabel('frequency shift from zero-bias cavity(kHz)')
plt.ylabel(r'Cross-correlation $S_{II_0}(\Delta f)$')


correlation_data = np.vstack((fAB, SAB))
np.savetxt(correlat_path+'correlation1mA.csv', correlation_data)



 # In[hist]
plt.plot(dat1[0], 10**(dat1[1]/20), label='zero bias, -30dB probe')
plt.plot(dat1[0], 10**(dat2[1]/20), label='1mA bias, -30dB probe')
plt.xlabel('frequency(GHz)')
plt.ylabel('lin mag')
plt.legend()

  # In[5] Plot density plot

#vlist = np.arange(0,len(filelist),1)
vlist = np.arange(0,11,1)
dat0 = path+files[0]
frequency = data_extract(data_folder, filelist[0])[0]

X,Y = np.meshgrid(frequency, vlist)

Z = np.zeros((len(vlist), len(frequency)))

for v in range(len(vlist)):
    dat = data_extract(Path=data_folder, file='100kO_vias{}V.s2p'.format(v))
    Z[v,:] = dat[1]

# In[5] Plot density plot
plt.pcolormesh(X,Y/100,Z)
plt.xlabel('Frequency(GHz)')
plt.ylabel('bias current(mA)')

cbar = plt.colorbar()
cbar.set_label('S21(dB)', rotation=270)
plt.tight_layout()

# In[5] Plot density plot

#vlist = np.arange(0,len(filelist),1)
vlist = np.arange(0,11,1)
dat0 = path+files[0]
frequency = data_extract(data_folder, filelist[0])[0]

X,Y = np.meshgrid(frequency, vlist)

Z = np.zeros((len(vlist), len(frequency)))

for v in range(len(vlist)):
    dat = data_extract(Path=data_folder, file='10kO_vias{}V.s2p'.format(v))
    Z[v,:] = dat[1]

# In[5] Plot density plot
plt.pcolormesh(X,Y/10,Z)
plt.xlabel('Frequency(GHz)')
plt.ylabel('bias current(mA)')

cbar = plt.colorbar()
cbar.set_label('S21(dB)', rotation=270)
plt.tight_layout()

# In[fit]

def lorentzian_asym(f, f0, df, Qe, Q, att, N):
    denom = sqrt(1 + 4*Q**2*((f-f0)/f0)**2)
    numerator = att*sqrt((1-Q/Qe)**2+4*Q**2*((f-f0)/f0-Q*df/(Qe*f0))**2)
    #denom = (1 + 4*Q**2*((f-f0)/f0)**2)
    #numerator = att*((1-Q/Qe)**2+4*Q**2*((f-f0)/f0-Q*df/(Qe*f0))**2)
    return N*numerator/denom

def inver(a, b):
    return a*b/(a-b)

def fitting_asym(file, guess, span1, span2):
    data = data_folder+file
    freq_raw = np.loadtxt(data, skiprows=23)[:,0]
    freq = freq_raw[np.where((freq_raw>guess[0]-span1) & (freq_raw <guess[0]+span2))[0]]
    s21_raw = np.loadtxt(data, skiprows=23)[:,3]
    s21 = s21_raw[np.where((freq_raw>guess[0]-span1) & (freq_raw <guess[0]+span2))[0]]
    #plt.plot(freq, 20*np.log10(s21), '.-')
    plt.plot(freq, 10**(s21/20), '.')

    popt, pcov = curve_fit(lorentzian_asym, freq, 10**(s21/20), p0=guess)
    plt.plot(freq, lorentzian_asym(freq, *popt), 'r-', label='fitted coupled Q = {}, \n fitted total Q ={}'.format(popt[2], popt[3])  )
    plt.plot(popt[0], lorentzian_asym(popt[0], *popt), 'go', markersize=15, label='reson freq = {}GHz'.format(popt[0]))
    plt.xlabel("frequency(GHz)")
    plt.ylabel("S21(lin mag)")
    plt.legend(loc='lower left')
    print(popt)
    print('reson freq = {}GHz,\n fitted coupled Q = {}, \n fitted total Q ={}, \n loss Q = {}'.format(popt[0], popt[2], popt[3], inver(popt[2], popt[3])))
    return(popt[0], popt[3], popt[2])
# In[actual fit]
fit_dat = fitting_asym(file=files, guess=[ 4.43280885e+00, -9.21322692e-06,  1.29053600e+04,  1.07487785e+04 ,2.25994734e+00,  7.64119378e-03], span1=0.005, span2=0.005)


# In[plot_fit]
filelist = os.listdir(data_folder)
reson_list = []

for files in filelist:
    data = data_extract(data_folder, files)
    fit_dat = fitting_asym(file=files, guess=[ 4.43280885e+00, -9.21322692e-06,  1.29053600e+04,  1.07487785e+04 ,2.25994734e+00,  7.64119378e-03], span1=0.005, span2=0.005)
    reson_list.append(fit_dat[0])

# In[fitted_reson_freq]
plt.figure(figsize=(10,5))  
bias_list=[0, 0.1, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0, 1/2, 0.1/2, 0.2/2, 0.3/2, 0.4/2, 0.5/2, 0.6/2, 0.7/2, 0.8/2, 0.9/2]  
plt.plot(bias_list, (np.array(reson_list)-reson_list[0])*1e6, '.')
plt.ylabel('fitted resonance frequency shift(kHz)')
plt.xlabel('bias current(mA)')

# In[plot for image]
plt.figure(figsize=(10,5))
fit_dat = fitting_asym(file=filelist[0], guess=[ 4.43280885e+00, -9.21322692e-06,  1.29053600e+04,  1.07487785e+04 ,2.25994734e+00,  7.64119378e-03], span1=0.005, span2=0.005)
'''
fitted data 

[ 4.43280926e+00 -3.79349294e-05  5.29322626e+04  2.06862253e+04
  2.15906112e+00  7.99556664e-03]
reson freq = 4.432809260881607GHz,
 fitted coupled Q = 52932.26261280987, 
 fitted total Q =20686.225251405696, 
 loss Q = 33956.69040518236
'''

# In[2022/11/30 : directory]
path = 'Z:/Mingi Kim/Resonators_0823~/1129_Leiden_mount/1130_biased resonator/'
data_folder = path + '-50dB/'

# In[meas-50dB]
'''
small SNR --> increase avg time to 3min (10sweep --> 30m.)

0V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 180, data_folder, '-50dB_10kO_vias0V.s2p')

# In[meas-50dB]
'''
1V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 180, data_folder, '-50dB_10kO_vias1V.s2p')

# In[meas-50dB]
'''
2V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 180, data_folder, '-50dB_10kO_vias2V.s2p')

# In[meas-50dB]
'''
3V
'''
S21trans(4.4322e9, 4.43325e9, 500, 100, 180, data_folder, '-50dB_10kO_vias3V.s2p')


# In[plot]
filelist = os.listdir(data_folder)
datlist = []

for files in filelist:
    data = data_extract(data_folder, files)
    plt.plot(data[0], data[1], '.')

plt.xlabel("frequency(GHz)")
plt.ylabel("S21(dB)")



