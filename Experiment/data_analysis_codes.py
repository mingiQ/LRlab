# In[1]:
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from scipy.optimize import curve_fit
from numpy import log10, pi, absolute, sqrt


# In[analysis codes]

def data_sorting():
    sweepfiles = os.listdir(wdir)

    sweeplist = [x for x in sweepfiles if "DAT" in x]
    fieldlist = []
    for files in sweeplist: fieldlist.append(float(files.split('_')[2].split('mT')[0]))

    #print(fieldlist)

    fieldlist, sweeplist = zipsort(fieldlist, sweeplist)
    
    return fieldlist, sweeplist

def dBmtoPhoton(power, freq):
    h = 6.62607015e-34
    n = 10**(power/10)/h*freq
    return n

def data_extract(Path=True, file=True, data_type=False, delimiter=','):    
    skip = 23
    
    if data_type == 'anri':
        skip = 23
    elif data_type == 'son':
        skip = 1
    elif data_type == 'rfsoc':
        skip = 2
    elif data_type == 'pna':
        skip = 8
        
    dat = np.loadtxt(Path+file, skiprows=skip)
    freq = dat[:,0]
    s21 = dat[:,3]
    phase = dat[:,4]
    return (freq, s21, phase)

'''
helpder functions for RFSoC phase calibration

'''
############################################################################################

def calculate_phase(d):
    [xi,xq] = d
    x = xi +1j*xq

    # Average to improve calibration.
    xavg = np.mean(x)

    # Calculate calibration phase.
    fi = np.remainder(np.angle(xavg,deg=True)+360,360)
    return [fi, np.abs(xavg), np.std(x)]


def phase_residuals(data,prediction):
    '''
    (data - pred + pi) mod 2pi ; returns the difference between data and pred
    '''
    r = np.remainder(data-prediction+180,360)-180
    return r
    
def phase_model(x, f0):
    '''
    x = [x0, x1] , 
    x0 : phase offset, x1 : ADC delay
    (x0 - 2pi*x1*f0) mod 2pi
    '''
    return np.remainder(x[0] - 360*x[1]*(f0), 360)

def phase_func(x, arg, freq):
    '''
    difference b/w data and phase drag due to the ADC delay = real signal
    '''
    resid = phase_residuals(arg, phase_model(x, freq))
    return resid

################################################################################################

def sweep_fit():
    '''
    fitting large amount of data - useful for field/power sweep
    todo : define parameters 
    
    '''
    
    sweepfiles = os.listdir(wdir_fieldsweep)

    sweeplist = [x for x in sweepfiles if "DAT" in x]
    fieldlist = []
    for files in sweeplist: fieldlist.append(float(files.split('_')[2].split('mT')[0]))

    #print(fieldlist)

    fieldlist, sweeplist = zipsort(fieldlist, sweeplist)

    flist = []
    qlist = []
    qinlist = []
    qexlist = []
    
    no = 23
    
    for files in sweeplist[:no]:
    
        data = data_extract(Path=wdir_fieldsweep, file=files, data_type='rfsoc')
        
        xdata = data[0][17:]
        ydata_raw = data[1][17:]
        ydata = 10*np.log10(ydata_raw)
        
        f0 = xdata[np.argmin(ydata_raw)]
        print(f0)
        
        title = files.split('_')[2]
        
        plt.figure(figsize=(10,5))
        p = fithanger_new_withQc(xdata, ydata, fitparams=[f0, 10000, 5000, 0, 50], domain=None, showfit=True, showstartfit=False,
                                 printresult=True, label=title+"mT", mark_data='bo', mark_fit='r-')
        
        
        Q_tot = (p[2]*p[1])/(p[2]+p[1])
        
        flist.append(p[0])
        qexlist.append(p[2])
        qinlist.append(p[1])
        
        qlist.append(Q_tot)
        
    fitsavedata = np.hstack((np.array([fieldlist[:no]]).T, np.array([flist]).T, np.array([qexlist]).T,np.array([qinlist]).T))
    np.savetxt(fname='//lrgrouphd2.physics.purdue.edu/usr/Mingi Kim/Resonators/230510_Leiden_Mount/fieldsweep_RFSoC/fits/BzvsQ_fit_results.DAT' , X=fitsavedata, delimiter= ',', header='Bz(mT), f0(MHz), Qex, Qint')
        
    plt.figure(figsize=(5,5))
    plt.plot(fieldlist[:no], flist, '.')
    plt.xlabel('Bz(mT)')
    plt.ylabel('Freq(MHz)')
    
    plt.figure(figsize=(10,5))
    
    ax1 = plt.subplot(1,2,1)
    ax1.plot(fieldlist[:no], qexlist, 'b.')
    ax1.set_xlabel('Bz(mT)')
    ax1.set_ylabel(r"$Q_{ex}$")
    
    ax2 = plt.subplot(1,2,2)
    ax2.plot(fieldlist[:no], qinlist, 'rx')
    ax2.set_xlabel('Bz(mT)')
    ax2.set_ylabel(r"$Q_{in}$")
    
    return fitsavedata


#################################################################################################
'''
all fitting codes will be moved to dsfit.py

'''

def lorentzian(f, f0, Qe, Q, N):
    return N-(Q/Qe)**2/(1 + 4*Q**2*((f-f0)/f0)**2)

def lorentzian_asym(f, f0, df, Qe, Q, att, N):
    
    '''
    Asymmetric Lorentzian fit : Journal of Applied Physics 111, 054510 (2012)

    '''
    denom = sqrt(1 + 4*Q**2*((f-f0)/f0)**2)
    numerator = att*sqrt((1-Q/Qe)**2+4*Q**2*((f-f0)/f0-Q*df/(Qe*f0))**2)
    #denom = (1 + 4*Q**2*((f-f0)/f0)**2)
    #numerator = att*((1-Q/Qe)**2+4*Q**2*((f-f0)/f0-Q*df/(Qe*f0))**2)
    return N*numerator/denom

def inver(a, b):
    return a*b/(a-b)

def fitting(path, file, guess, span1, span2):
    data = path+file
    freq_raw = np.loadtxt(data, skiprows=23)[:,0]
    freq = freq_raw[np.where((freq_raw>guess[0]-span1) & (freq_raw <guess[0]+span2))[0]]
    s21_raw = np.loadtxt(data, skiprows=23)[:,3]
    s21 = s21_raw[np.where((freq_raw>guess[0]-span1) & (freq_raw <guess[0]+span2))[0]]
    #plt.plot(freq, 20*np.log10(s21), '.-')
    plt.plot(freq, 10**(s21/20), '.')

    popt, pcov = curve_fit(lorentzian, freq, 10**(s21/20), p0=guess) 
    
    plt.plot(freq, lorentzian(freq, *popt),'-', label='reson freq = {}GHz,\n fitted coupled Q = {}, \n fitted total Q ={}'.format(popt[0], popt[1], popt[2]) )
    plt.legend(loc='upper right')
    
    print('reson freq = {}GHz,\n fitted coupled Q = {}, \n fitted total Q ={}, \n loss Q = {}'.format(popt[0], popt[1], popt[2], inver(popt[1], popt[2])))


def fitting_asym(data_folder, file, guess, span1, span2):
    
    '''
    data_folder : working directory
    file: filename
    guess : [freq, df, external Q, total Q, scale1, scale2]
    span1 : fit region left
    span2 : fit region right
    
    '''
    data = data_folder+file
    freq_raw = np.loadtxt(data, skiprows=23)[:,0]
    freq = freq_raw[np.where((freq_raw>guess[0]-span1) & (freq_raw <guess[0]+span2))[0]]
    s21_raw = np.loadtxt(data, skiprows=23)[:,3]
    s21 = s21_raw[np.where((freq_raw>guess[0]-span1) & (freq_raw <guess[0]+span2))[0]]
    #plt.plot(freq, 20*np.log10(s21), '.-')
    plt.plot(freq, 10**(s21/20), '.')

    popt, pcov = curve_fit(lorentzian_asym, freq, 10**(s21/20), p0=guess, maxfev=5000)
    plt.plot(freq, lorentzian_asym(freq, *popt), 'r-', label='fitted coupled Q = {}, \n fitted total Q ={}'.format(popt[2], popt[3])  )
    plt.plot(popt[0], lorentzian_asym(popt[0], *popt), 'go', markersize=15, label='reson freq = {}GHz'.format(popt[0]))
    plt.xlabel("frequency(GHz)")
    plt.ylabel("S21(lin mag)")
    plt.legend(loc='lower left')
    print(popt)
    print('reson freq = {}GHz,\n fitted coupled Q = {}, \n fitted total Q ={}, \n loss Q = {}'.format(popt[0], popt[2], popt[3], inver(popt[2], popt[3])))
    return(popt[0], popt[3], popt[2])

##############################################################################################

def cross_correlate(dat1, dat2, step, scan_range, correlat_path, filename):

    sum_range = np.arange(-scan_range,scan_range+1,1)
    freq_renorm = np.linspace(dat1[0][0], dat1[0][-1], 4001)
    # step = 2.625e-7  # 0.26kHz
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

def square_correlate(dat1, dat2, scan_range, correlat_path, filename):

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
    

def densityplot_2d(x_data, y_data, z_data, colormap, label_x, label_y, label_z):
     X,Y = np.meshgrid(x_data, y_data)
     
     Z = np.zeros((len(y_data), len(x_data)))
 
     for y in range(len(y_data)):
     
         Z[y,:] = z_data[y]
 
     plt.pcolormesh(X,Y,Z, cmap=colormap)
     plt.xlabel('Frequency(GHz)')
     plt.ylabel('power(dBm)')
     plt.colorbar(label="S21(dB)")
    

# In[correlation fucntion between spectra]
# =============================================================================
# 
# dat1 = data_extract(data_folder, filelist[22])  # reference signal
# dat2 = data_extract(data_folder, filelist[23])  # test signal
# freq_renorm = np.linspace(dat1[0][0], dat1[0][-1], 4001)
# 
# 
# slices = len(dat1[0])
# 
# '''
# say, integration region to be 4.4325GHz to 4.4331GHz (linewidth ~ 200kHz) larger than cavity linewidth
# 
# I matched frequency step with 262.5Hz / total dataset starts from 4.4322GHz to 4.43325GHz, 4001 slices, (index runs from 0 to 4000)
# '''
# step = 2.625e-7
# # =============================================================================
# # 
# # df= []
# # 
# # for i in range(slices-1):
# #     df.append(freq_renorm[i+1]-freq_renorm[i])
# # 
# # =============================================================================
# 
# '''
# cross-correlation function at f=-1000*step
# 
# '''
# fAB=[]
# SAB=[]
# 
# fAB.append(-1000*step)
# #plt.plot(freq_renorm, 10**(dat1[1]/20))  # ref
# #plt.plot(freq_renorm+1000*step, 10**(dat2[1]/20))
# #plt.plot(freq_renorm-1000*step, 10**(dat2[1]/20)) # input
# 
# # SB[f+DF] = freq_renorm-1000*step, 10**(dat2[1]/20)
# '''
# SA --> summation ranges over 1000 ~ 3000
# SB --> summation ranges over 2000 ~ 4000
# 
# '''
# plt.plot(freq_renorm[1000:3000], (10**(dat1[1]/20))[1000:3000])  # ref
# #plt.plot(freq_renorm+1000*step, 10**(dat2[1]/20))
# plt.plot((freq_renorm-1000*step)[2000:4000], (10**(dat2[1]/20))[2000:4000]) # input
# 
# s=0
# for i in range(len(freq_renorm[1000:3000])):
#    s+= (10**(dat1[1]/20))[1000:3000][i]*(10**(dat2[1]/20))[2000:4000][i]*step
#    
# SAB.append(s)
# 
# # In[cross-correlation function]
# sum_range = np.arange(-1000,1001,1)
# step = 2.625e-7
# fAB=[]
# SAB=[]
# 
# for i in range(len(sum_range)):
#     order = sum_range[i]
#     fAB.append(order*step)
#     s = 0
#     for ii in range(len(freq_renorm[1000:3000])):
#         s+= (10**(dat1[1]/20))[1000:3000][ii]*(10**(dat2[1]/20))[1000-order:3000-order][ii]*step
#         
#     SAB.append(s)
# 
# # In[cross corrleation plot]
# correlat_path = 'C:/Users/20141/OneDrive - purdue.edu/PhD_research/resonators/1129_Leiden_mount/1130_biased resonator/correlaton function/'
# 
# plt.plot(np.array(fAB)*1e6, SAB, '.')
# plt.xlabel('frequency shift from zero-bias cavity(kHz)')
# plt.ylabel(r'Cross-correlation $S_{II_0}(\Delta f)$')
# 
# 
# correlation_data = np.vstack((fAB, SAB))
# np.savetxt(correlat_path+'correlation1mA.csv', correlation_data)
# 
# 
# 
#  # In[hist]
# plt.plot(dat1[0], 10**(dat1[1]/20), label='zero bias, -30dB probe')
# plt.plot(dat1[0], 10**(dat2[1]/20), label='1mA bias, -30dB probe')
# plt.xlabel('frequency(GHz)')
# plt.ylabel('lin mag')
# plt.legend()
# 
#   # In[5] Plot density plot
# 
# #vlist = np.arange(0,len(filelist),1)
# vlist = np.arange(0,11,1)
# dat0 = path+files[0]
# frequency = data_extract(data_folder, filelist[0])[0]
# 
# X,Y = np.meshgrid(frequency, vlist)
# 
# Z = np.zeros((len(vlist), len(frequency)))
# 
# for v in range(len(vlist)):
#     dat = data_extract(Path=data_folder, file='100kO_vias{}V.s2p'.format(v))
#     Z[v,:] = dat[1]
# 
# # In[5] Plot density plot
# plt.pcolormesh(X,Y/100,Z)
# plt.xlabel('Frequency(GHz)')
# plt.ylabel('bias current(mA)')
# 
# cbar = plt.colorbar()
# cbar.set_label('S21(dB)', rotation=270)
# plt.tight_layout()
# 
# # In[5] Plot density plot
# 
# #vlist = np.arange(0,len(filelist),1)
# vlist = np.arange(0,11,1)
# dat0 = path+files[0]
# frequency = data_extract(data_folder, filelist[0])[0]
# 
# X,Y = np.meshgrid(frequency, vlist)
# 
# Z = np.zeros((len(vlist), len(frequency)))
# 
# for v in range(len(vlist)):
#     dat = data_extract(Path=data_folder, file='10kO_vias{}V.s2p'.format(v))
#     Z[v,:] = dat[1]
# 
# # In[5] Plot density plot
# plt.pcolormesh(X,Y/10,Z)
# plt.xlabel('Frequency(GHz)')
# plt.ylabel('bias current(mA)')
# 
# cbar = plt.colorbar()
# cbar.set_label('S21(dB)', rotation=270)
# plt.tight_layout()
# 
# 
# 
# # In[visualization codes]
# 
# def densityplot_2d(x_data, y_data, z_data, colormap, label_x, label_y, label_z):
#     X,Y = np.meshgrid(x_data, y_data)
#     
#     Z = np.zeros((len(y_data), len(x_data)))
# 
#     for y in range(len(y_data)):
#     
#         Z[y,:] = z_data[y]
# 
#     plt.pcolormesh(X,Y,Z, cmap='inferno')
#     plt.xlabel('Frequency(GHz)')
#     plt.ylabel('power(dBm)')
#     plt.colorbar(label="S21(dB)")
# =============================================================================
