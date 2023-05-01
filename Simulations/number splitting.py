# In[0]: pkgs

import numpy as np
from numpy import tan, pi, sqrt
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import minimize, curve_fit
from scipy.constants import e,h,hbar,k,G,epsilon_0

phi_0 = hbar/2./e
from scipy import signal
import matplotlib as mpl
import pickle

from qutip import *
qutip.settings.has_mkl=False
#%matplotlib notebook

#from lmfit import Model, Parameters

# In[1]

def spectroscopy_fit(fd_list, eps_q, eps_c, gamma, alpha, beta):
    
    Nq = 5            # Dimension of qubit.
    Nc = 7            # Dimension of resonator. 
    fq = 4.75e9      # Frequency of qubit.
    fc = 3.691e9      # Frequency of resoantor.
    Anh = 245e6       # Qubit anharmonicity
    chi = 9e6         # Dispersive shift 
    kappa = 0.5e6     # Resoantor linewidth.
    


# Fitting parameters.

#     alpha     : Scale factor.
#     beta      : Offset.
#     gamma     : Qubit natural linewidth.
#     eps_q     : Qubit Drive amplitude.
#     eps_c     : Qubit Drive amplitude.
    
    a = tensor(destroy(Nq),qeye(Nc)) 
    b = tensor(qeye(Nq),destroy(Nc)) 
    
    num_a = a.dag()*a
    num_b = b.dag()*b
    
    fq = 2*pi*fq
    fc = 2*pi*fc
    
    delta = np.abs(fc-fq)
    g = np.sqrt(chi*delta)
    
    Anh = 2*pi*Anh
    chi = 2*pi*chi
    eps_q = 2*pi*eps_q
    eps_c = 2*pi*eps_c
    gamma = 2*pi*gamma
    kappa = 2*pi*kappa
    
    fp = fc # Probe frequency in resonator resonance

    r=[]
    
    H0 = (fq)*num_a + (fc)*num_b - 0.5*Anh*a.dag()*a.dag()*a*a - chi*num_a*num_b
    
    for fd in fd_list:                             # qubit pump frequencies.
        
        fd = 2*pi*fd
        H = H0 
        
        H -= fd*num_a #+ fp*num_a
        H -= fp*num_b #+ fd*num_b
        H += +eps_q*(b.dag()+b) + eps_c*(b.dag()+b)
        H += +eps_c*g/delta*(a.dag()+a) + eps_q*g/delta*(a.dag()+a)
                
        c_ops = []
        c_ops.append(np.sqrt(gamma)*a)
        c_ops.append(np.sqrt(kappa)*b)
        
        rho_ss = steadystate(H, c_ops)
        r.append(expect(num_a,rho_ss))
        
    return alpha*np.array(r) + beta  

# In[2]:
# Initial guess.
alpha = .1
beta  = 1.2e-5

fd_list = np.linspace(4.7,4.77,101)*1e9
eps_c = 100e3
eps_q = 100e3
gamma = 225e3

sim = spectroscopy_fit(fd_list,eps_q, eps_c, gamma, alpha, beta)

#plt.figure()
plt.plot(fd_list, sim, '.-',label='Model')
plt.xlabel("$f$ (GHz)")
plt.ylabel("Spec signal")


# In[3]:
    
#