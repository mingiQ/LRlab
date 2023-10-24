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

# In[Andreev level visualization]

Delta = 20.8e9 # GHz
tau = 0.96   # balistic junction

phi_list = np.linspace(0, 2, 100)
plt.plot(phi_list, [E_A(Delta, tau, ph*pi)/1e9 for ph in phi_list], label=r'$\left| e \right\rangle$')
plt.plot(phi_list, [-E_A(Delta, tau, ph*pi)/1e9 for ph in phi_list], label=r'$\left| g \right\rangle$')
plt.plot(phi_list, [E_A(Delta, 1, ph*pi)/1e9 for ph in phi_list], 'k--')
plt.plot(phi_list, [-E_A(Delta, 1, ph*pi)/1e9 for ph in phi_list], 'k--')
plt.axhline(y = 0, color = 'r', linestyle = '-', label=r'$\left| o \right\rangle$')
plt.title(r"$\tau$={}, $2E_A$={}GHz @ $\varphi=\pi$".format(tau, np.round(2*E_A(Delta, tau,pi)/1e9,1) ))
plt.legend(loc="lower center")
plt.xlabel(r"$\phi/\pi$")
plt.ylabel("E/h(GHz)")

# In[Andreev current visualization]

Delta = 20.8e9 # GHz
tau = 0.98   # balistic junction

phi_list = np.linspace(-1, 3, 100)
plt.plot(phi_list, [0.5*I_A(Delta, tau, ph*pi)[0]*1e9 for ph in phi_list], label=r'$\left| e \right\rangle$')
plt.plot(phi_list, [0.5*I_A(Delta, tau, ph*pi)[1]*1e9 for ph in phi_list], label=r'$\left| g \right\rangle$')

plt.plot(phi_list, [0.5*I_A(Delta, 1, ph*pi)[0]*1e9 for ph in phi_list], 'k--')
#plt.plot(phi_list, [-E_A(Delta, 1, ph*pi)/1e9 for ph in phi_list], 'k--')
plt.axhline(y = 0, color = 'r', linestyle = '-', label=r'$\left| o \right\rangle$')
plt.title(r"$\tau$={}, $2E_A$={}GHz @ $\varphi=\pi$".format(tau, np.round(2*E_A(Delta, tau,pi)/1e9,1) ))
plt.legend(loc="lower center")
plt.xlabel(r"$\phi/\pi$")
plt.ylabel(r"$I_A$(nA)")

# In[Jaynes-Cummings approx.]
M = 10e-12
Zr = 50
fr = 6e9

E = eigen_ABS_JC(1, tau, np.pi, M, Zr, fr, Delta)
nph = 0


ph = np.linspace(0.9, 1.1, 100)
plt.figure(figsize=(8,8))
plt.title("ABS-cavity photon eigenenergy soultion\n from Jaynes-Cummings Hamiltonian")

plt.plot(ph, [eigen_ABS_JC(nph, tau, x*np.pi, M, Zr, fr, Delta) for x in ph ])
plt.plot(ph, [eigen_ABS_JC(nph+1, tau, x*np.pi, M, Zr, fr, Delta) for x in ph ])
plt.annotate(r'photon mode : {}'.format(nph)+'\n'+
             r'$\tau$ : {}'.format(tau)+' \n'+
             ' M : {} pH '.format(M*1e12)+'\n'+
             r' $Z_R$= {} $\Omega$'.format(Zr) +'\n'+
             r' $f_R$ = {}GHz'.format(fr/1e9), xy = (1, 0.9), xycoords='axes fraction')

plt.xlabel(r"reduced phase $\varphi/\pi$")
plt.ylabel(r"transition energy ($E_t/f_r$)")
# In[]
ph = np.linspace(0.97, 1.03, 100)
plt.figure(figsize=(8,8))
plt.title("ABS-cavity photon eigenenergy soultion\n from Jaynes-Cummings Hamiltonian")

plt.plot(ph, [eigen_ABS_JC(nph, tau, x*np.pi, M, Zr, fr, Delta) for x in ph ])
#plt.plot(ph, [eigen_ABS_JC(nph+1, tau, x*np.pi, M, Zr, fr, 45e9) for x in ph ])
plt.annotate(r'photon mode : {}'.format(nph)+'\n'+
             r'$\tau$ : {}'.format(tau)+' \n'+
             ' M : {} pH '.format(M*1e12)+'\n'+
             r' $Z_R$= {} $\Omega$'.format(Zr) +'\n'+
             r' $f_R$ = {}GHz'.format(fr/1e9), xy = (1, 0.9), xycoords='axes fraction')

plt.xlabel(r"reduced phase $\varphi/\pi$")
plt.ylabel(r"transition energy ($E_t/f_r$)")

# In[]

plt.plot(ph, [eigen_ABS_JC(nph, tau, x*np.pi, M, Zr, fr, 45e9)[0]-eigen_ABS_JC(nph, tau, x*np.pi, M, Zr, fr, 45e9)[1] for x in ph ])

