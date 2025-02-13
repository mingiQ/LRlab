# In[1]:
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from scipy.optimize import curve_fit
from numpy import log10, pi, absolute, sqrt
from scipy.constants import h,hbar, e, Boltzmann

# In[analysis codes]

def data_extract(Path, file, data_type):    
    skip = 23
    
    if data_type == 'anri':
        skip = 23
    elif data_type == 'son':
        skip = 1
    elif data_type == 'rfsoc':
        skip = 2
    elif data_type == 'pna':
        skip = 8
    elif data_type == 'cmt':
        skip = 5
        
    dat = np.loadtxt(Path+file, skiprows=skip)
    freq = dat[:,0]
    s21 = dat[:,3]
    phase = dat[:,4]
    return (freq, s21, phase)

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


# =============================================================================
# from PyQt5 import QtWidgets, QtCore
# from pyqtgraph import PlotWidget, plot
# import pyqtgraph as pg
# import sys  # We need sys so that we can pass argv to QApplication
# import os
# 
# class PlotWindow(QtWidgets.QMainWindow):
# 
#     def __init__(self, xs, ys, *args, **kwargs):
#         super(PlotWindow, self).__init__(*args, **kwargs)
# 
#         self.graphWidget = pg.PlotWidget()
#         self.setCentralWidget(self.graphWidget)
# 
#         #hour = [1,2,3,4,5,6,7,8,9,10]
#         #temperature_1 = [30,32,34,32,33,31,29,32,35,45]
#         #temperature_2 = [50,35,44,22,38,32,27,38,32,44]
# 
#         #Add Background colour to white
#         self.graphWidget.setBackground('w')
#         # Add Title
#         self.graphWidget.setTitle("Your Title Here", color="b", size="30pt")
#         # Add Axis Labels
#         styles = {"color": "#f00", "font-size": "20px"}
#         self.graphWidget.setLabel("left", "Transmission(ADC levels)", **styles)
#         self.graphWidget.setLabel("bottom", "Freq(MHz)", **styles)
#         #Add legend
#         self.graphWidget.addLegend()
#         #Add grid
#         self.graphWidget.showGrid(x=True, y=True)
#         #Set Range
#         #self.graphWidget.setXRange(0, 10, padding=0)
#         #self.graphWidget.setYRange(20, 55, padding=0)
#         
#         for i in range(len(xs)):
#             self.plot(xs[i], ys[i], "{}th".format(i), 'r')
# 
#         #self.plot(hour, temperature_1, "Sensor1", 'r')
#         #self.plot(hour, temperature_2, "Sensor2", 'b')
# 
#     def plot(self, x, y, plotname, color):
#         pen = pg.mkPen(color=color)
#         self.graphWidget.plot(x, y, name=plotname, pen=pen, symbol='+', symbolSize=10, symbolBrush=(color))    
# =============================================================================

def densityplot_2d(x_data, y_data, z_data, colormap, label_x, label_y, label_z):
     X,Y = np.meshgrid(x_data, y_data)
     
     Z = np.zeros((len(y_data), len(x_data)))
 
     for y in range(len(y_data)):
     
         Z[y,:] = z_data[y]
 
     plt.pcolormesh(X,Y,Z, cmap='inferno')
     plt.xlabel('Frequency(GHz)')
     plt.ylabel('power(dBm)')
     plt.colorbar(label="S21(dB)")

# =====================================================================================================
'''
Andreev bound state analysis

'''


def z(M, Zr, fr):
    Rq = h/(4*e**2)
    omega_r = 2*pi*fr
    return pi*M**2*omega_r**2/(Zr*Rq)

def E_A(Delta_sc, tau, phi):
    E = Delta_sc*np.sqrt(1-tau*np.sin(phi/2)**2)
    return E

def I_A(Delta_sc, tau, phi):
    I = h*Delta_sc/(4*(h/(2*e)))*tau*np.sin(phi)/np.sqrt(1-tau*np.sin(phi/2)**2)
    return -I, I

def gc(tau, phi, M, Zr, fr, Delta_sc):
    E_a = E_A(Delta_sc, tau, phi)
    E_a_pi = E_A(Delta_sc, tau, np.pi)
    return np.sqrt(z(M, Zr, fr))*E_a_pi/2*(Delta_sc/E_a - E_a/Delta_sc)

def eigen_ABS_JC(n_ph, tau, phi, M, Zr, fr, Delta_sc):
    d = 2*E_A(Delta_sc, tau, phi) - fr
    Ep = (E_A(Delta_sc, tau, phi) + (n_ph+1/2)*fr + np.sqrt(gc(tau, phi, M, Zr, fr, Delta_sc)**2+d**2/4))/(fr)
    Em = (E_A(Delta_sc, tau, phi) + (n_ph+1/2)*fr - np.sqrt(gc(tau, phi, M, Zr, fr, Delta_sc)**2+d**2/4))/(fr)
    return Ep, Em

gamma = 0.00026 * 2*pi
kappa = 0.0013 * 2*pi
fr=9
w = fr * 2*pi
delta = 0.4
tau = 0.9999

M = 29e-12
Delta_sc = 320e9
def transmission_rfsquid_hanger(gamma, kappa, fr, w, delta, tau, Zr):
    f_lists = np.linspace(w/(2*pi)-delta/2, w/(2*pi)+delta, 800)
    phi_lists = np.linspace(0.97, 1.03, 800)
    
    t_array = np.zeros((len(f_lists), len(phi_lists)))
    for i, f in enumerate(f_lists):
        t_calc = []
        w_d = f * 2*pi
        for j, ph in enumerate(phi_lists):
            
            phi = ph*pi
            
            wq = 2*E_A(Delta_sc=Delta_sc, tau=tau, phi=phi)/1e9 * 2* pi
            g = gc(tau=tau, phi=phi, M=M, Zr=Zr, fr=fr*1e9, Delta_sc=Delta_sc)/1e9 * 2*pi
            
            #print(f"phase = {ph}, Andreev level: {wq/(2*pi)} GHz, coupling {g/(2*pi)} GHz")
            
            detuning_r = w_d - w + 1j*kappa/2
            detuning_q = w_d - wq + 1j*gamma
            
            t = 1-np.pi/2*1j*kappa/2/(detuning_r - g**2/detuning_q)
            #print(abs(t))
            t_calc.append(abs(t))
        t_array[i] = np.array(t_calc)
    
    #t_arr = t_array[::-1]
    plt.rcParams['font.size'] = '16'
    
    plt.title(fr'$\tau=${tau}, M={M*1e12}pH')
    
    plt.pcolor(phi_lists, (f_lists-fr)*1e3, t_array, cmap='RdBu')
    #plt.ylim(0.9,1.1)
    plt.ylabel('detuning (MHz)')
    plt.xlabel(r'phase bias ($\varphi/\pi$)')
    plt.colorbar(label='Normalized Transmission')  # Add colorbar with label
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
