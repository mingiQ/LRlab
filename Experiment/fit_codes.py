# In[0] : necessary pkgs

'''
fitting pkg ; ref: slab/slab/dsfit.py
https://github.com/SchusterLab/slab/blob/master/slab/dsfit.py
author : David Schuster
commented by : Mingi

'''
import numpy as np
import math as math
import matplotlib.pyplot as plt
import scipy
import scipy.fftpack
import cmath
import numpy

from scipy import optimize
from scipy.constants import h,hbar, e, Boltzmann
from numpy import pi

# In[1] : sorting data

# def set_fit_plotting(pkg='matplotlib'):
#     global plt
#     plt={'guiqwt':plt1,'matplotlib':plt2}[pkg]
    
def argselectdomain(xdata,domain):
    ind=np.searchsorted(xdata,domain)    
    # if domain is [f_start, f_stop], searchsorted returns the index of xdata where xdata[index-1] < f_i <= xdata[index] (left) 
    return (ind[0],ind[1])


def selectdomain(xdata,ydata,domain):
    '''
    choose certain region in xdata and correponding ydata
    xdata = array
    ydata = array
    domain = float or array e.g. [2.3e9, 4.5e9] : select x data where x[i] is larger than 2.3e9 and smaller than 4.5e9
    '''
    ind=np.searchsorted(xdata,domain)
    return xdata[ind[0]:ind[1]],ydata[ind[0]:ind[1]]

def zipsort(xdata,ydata):
    inds=np.argsort(xdata)
    return np.take(xdata,inds),np.take(ydata,inds,axis=0)
    
# In[2]: models

'''
Polynoimals
'''

def linear(p,x):
    return p[0]+p[1]*(x)

def poly(p, x):
    return p[1]*(x-p[-1])+p[2]*(x-p[-1])**2

def polynomial2(p,x):
    return p[0]+p[1]*(x-p[-1])+p[2]*(x-p[-1])**2

def polynomial(p, x):
    return p[0] + p[1]*(x - p[-1]) + p[2]*(x - p[-1])**2 + p[3]*(x - p[-1])**3 + p[4]*(x - p[-1])**4 + \
           p[5]*(x - p[-1])**5 + p[6]*(x - p[-1])**6 + p[7]*(x - p[-1])**7 + p[8]*(x - p[-1])**8 + p[9]*(x - p[-1])**9


'''
exponentials
'''

def expfunc(p,x):
    """p[0]+p[1]*exp(-(x-p[2])/p[3])"""
    return p[0]+p[1]*math.e**(-(x-p[2])/p[3])

def doubleexpfunc(p, x):
    return (p[0]+p[1]*math.e**(-(x-p[2])/p[3])+ p[4]*math.e**(-(x-p[2])/10000))
    # Only one offset

def gaussfunc_nooffset(p, x):
    """p[0] exp(- (x-p[1])**2/p[2]**2/2)"""
    return p[0]*math.e**(-1./2.*(x-p[1])**2/p[2]**2)

def gaussfunc(p, x):
    """p[0]+p[1] exp(- (x-p[2])**2/p[3]**2/2)"""
    return p[0]+p[1]*math.e**(-1./2.*(x-p[2])**2/p[3]**2)
    


'''
Transmissions
'''

def hangerfunc_old(p,x):
    """p=[f0,Q,S21Min,Tmax]
       (4*(x-p[0])**2/p[0]**2 * p[1]**2+p[2]**2/(1+4*(x-p[0])**2/p[0]**2 * p[1]**2))*p[3]  linmag
    """
    return ((4.*((x-p[0])* p[1]/p[0])**2. +p[2]**2.)/(1.+4.*((x-p[0])* p[1]/p[0])**2.))*p[3]

def hangerfunc(p,x):
    """p=[f0,Qi,Qc,df,scale]  linmag
        ! different from [An analysis method for asymmetric resonator transmission applied to superconducting devices]
        Journal of Applied Physics 111, 054510 (2012) https://doi.org/10.1063/1.3692073
    """
    #print p    
    f0,Qi,Qc,df,scale = p
    a=(x-(f0+df))/(f0+df)
    b=2*df/f0
    Q0=1./(1./Qi+1./Qc)
    return scale * (-2. * Q0 * Qc + Qc ** 2. + Q0 ** 2. * (1. + Qc ** 2. * (2. * a + b) ** 2.)) / (
    Qc ** 2 * (1. + 4. * Q0 ** 2. * a ** 2.))

def hangerfunc_new(p,x):
    """p=[f0,Qi,Qc,df,scale]  logmag """
    #print p    
    f0,Qi,df,scale = p
    a=(x-(f0+df))/(f0+df)
    b=2*df/f0
    Qc=4000.
    y = 10 * np.log10(scale * (Qc ** 2. + Qi ** 2. * Qc ** 2. * (2. * a + b) ** 2.) / (
    (Qc + Qi) ** 2 + 4. * Qi ** 2. * Qc ** 2. * a ** 2.))
    return y
    
def hangerfunc_new_withQc(p,x):
    """p=[f0,Qi,Qc,df,scale]   logmag -- similar  to JAP111, 054510 but not exactly the same"""
    #print p    
    f0,Qi,Qc,df,scale = p
    a=(x-(f0+df))/(f0+df)
    b=2*df/f0
    y = 10 * np.log10(scale * (Qc ** 2. + Qi ** 2. * Qc ** 2. * (2. * a + b) ** 2.) / (
    (Qc + Qi) ** 2 + 4. * Qi ** 2. * Qc ** 2. * a ** 2.))
    return y

def hangerfunc_circle(p,x):
    '''
    Parameters
    ----------
    p : []
    x : frequency
    Returns
    -------
    None.

    '''

def hangerfunctilt(p,x):
    """Ge Editing  p=[f0,Qi,Qc,df,scale,slope, offset] linmag"""
    f0, Qi, Qc, df, slope, offset = p
    a  = (x-(f0+df))/(f0+df)
    b  = 2*df/f0
    Q0 = 1./(1./Qi+1./Qc)
    y  = math.exp(slope*x+offset)
    #y=[math.exp(slope*i+offset) for i in x]
#    return slope*x+offset+scale*(-2.*Q0*Qc + Qc**2. + Q0**2.*(1. + Qc**2.*(2.*a + b)**2.))/(Qc**2*(1. + 4.*Q0**2.*a**2.))
    return y * (-2. * Q0 * Qc + Qc ** 2. + Q0 ** 2. * (1. + Qc ** 2. * (2. * a + b) ** 2.)) / (
    Qc ** 2 * (1. + 4. * Q0 ** 2. * a ** 2.))

def lorentzian_asym(p, f):
    
    '''
    Asymmetric Lorentzian fit : Journal of Applied Physics 111, 054510 (2012)

    '''
    f0, df, Qe, Q, scale = p
    denom = np.sqrt(1 + 4*Q**2*((f-f0)/f0)**2)
    numerator = scale * np.sqrt((1-Q/Qe)**2+4*Q**2*((f-f0)/f0-Q*df/(Qe*f0))**2)
    return numerator/denom


'''
ETC
'''

def harmfunc(p, x):
    """p[0]+p[1]/(1+(x-p[2])**2/p[3]**2)"""
    return p[0]+p[1]/((x**2-(p[2])**2)**2 + (p[3]**2)*x**2)**(0.5)

def lorfunc(p, x):
    """p[0]+p[1]/(1+(x-p[2])**2/p[3]**2)"""
    return p[0]+p[1]/(1+(x-p[2])**2/p[3]**2)


    
def pulse_errfunc(p,x):
    """p[0]+p[1]*exp(-(x-p[2])/p[3])"""
    return p[0]+0.5*(1-((1-p[1])**x))

def decaysin(p,x):
    """p[0]*np.sin(2.*pi*p[1]*x+p[2]*pi/180.)*np.e**(-1.*(x-p[5])/p[3])+p[4]"""
    return p[0]*np.sin(2.*np.pi*p[1]*x+p[2]*np.pi/180.)*np.e**(-1.*(x-p[5])/p[3])+p[4]

def SNT_func(p, v):
    Tn, GB, T, voff = p
    qe, kb = (1.6e-19, 1.38e-23)
    return GB * kb * (Tn + .5 * (qe * (v - voff) / kb) / np.tanh(qe * (v - voff) / (2 * kb * T)))



def decayrabifunc1(p, x):
    return (p[0]+ p[1]*(math.e**(-(x-p[2])/p[3]))*((np.cos(np.pi*p[4]*(x-p[5])))**2))
    # Only one offset

def rabiwidth(p, x):
    return np.sqrt(p[1]**2*(x**2) + p[0]**2)

def rabiwidth2(p, x):
    return np.sqrt((p[1]-x)**2 + p[0]**2)

def rabisatfunc(p, x):
    return p[0] + x**2/(2*x**2 + 1/p[1] )

def dispersiveshift(p, x):
    """p[0]+p[1]/(1+(x-p[2])**2/p[3]**2)"""
    return p[0]+p[1]/(1+(x-p[2])**2/p[3]**2) + p[4]/(1+(x-p[5])**2/p[6]**2)

'''
ABS spectrum in short junction
'''

def z(M, Zr, fr):
    Rq = h/(4*e**2)
    omega_r = 2*pi*fr
    return pi*M**2*omega_r**2/(Zr*Rq)

def E_A(p, phi):
    Delta_sc, tau = p
    E = Delta_sc*np.sqrt(1-tau*np.sin(phi/2)**2)
    return E

def I_A(p, phi):
    Delta_sc, tau = p
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



## TLS fit ref: arxiv 1010.6063v2 eq 5, 6

def insertion_loss(Qint, Qex, Pinc, n):
    '''
    calculates insertion loss of coupling capacitor (or inductor)
    and circulating power inside the resonator. 
    Pinc: incident power right before the resonator. Probe power / line attenuation --> lin mag!!
    n: nth harmonics of the resonator 

    output: 
    IL (insertion loss) in lin mag
    Pcirc : circulating power
    '''
    Qtot = 1/(1/Qint+1/Qex)  # loaded Q
    IL = 1-Qtot/Qint  # insertion loss in lin mag
    Pcirc = Qtot/(n*np.pi) * (IL * Pinc)
    return IL, Pcirc

def tlsfit_model(p, x):
    '''
    we will fit inverse TLS function
    ref: arxiv 1010.6063v2 eq 5, 6
    
    x : converted circulating power e.g. power_circ = np.array([insertion_loss(Qint[n], Qex[n], power_in_W[n], 1)[1] for n in range(len(power_in_W))])
    '''
    h = 6.62607015e-34 
    kB = 1.380649e-23
    Q0, losstan, beta, Pc, f0, T = p
    inverseQtls = losstan * np.tanh(h*f0/(2*kB*T))/np.sqrt(1+(x/Pc)**beta)
    y = inverseQtls + 1/Q0
    return y

# In[3]: extract parameters

def print_cavity_Q(fit):
    print(fit[2]/2/fit[3])
    return fit[2]/2/fit[3]


def hangerqs_old(fitparams):
    """Converts fitparams into Qi and Qc"""
    return abs(fitparams[1]/fitparams[2]), abs(fitparams[1])/(1-abs(fitparams[2]))


# In[4]: fit
    
"""Wraplter around scipy.optimize.leastsq"""


def fitgeneral(xdata, ydata, fitfunc, fitparams, domain=None, showfit=False, showstartfit=False, showdata=True,
               label="", mark_data='bo', mark_fit='r-'):
    """Uses optimize.leastsq to fit xdata ,ydata using fitfunc and adjusting fit params"""
    '''
    xdata : 1d array
    ydata : 1d array
    fitfunc : above model functions in In[2]
    fitparams : inital guess for the fitting parameters p_init = [p0,p1,...] ~ 1d array
    domain : [x_start, x_stop]
    
    '''

    # sort data
    order = np.argsort(xdata)  
    xdata = xdata[order]
    ydata = ydata[order]

    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
        
#    print 'minimum', np.min(fitdatay)
#    ymin=np.min(fitdatay)
    errfunc = lambda p, x, y: (fitfunc(p,x) - y) #there shouldn't be **2 # Distance to the target function
    
    startparams=fitparams # Initial guess for the parameters
    
    bestfitparams, success = optimize.leastsq(errfunc, startparams[:], args=(fitdatax,fitdatay))
    
    if showfit:
        if showdata:
            plt.plot(fitdatax,fitdatay,mark_data,label=label+" data")
        if showstartfit:
            plt.plot(fitdatax,fitfunc(startparams,fitdatax),label=label+" startfit")
        plt.plot(fitdatax,fitfunc(bestfitparams,fitdatax),mark_fit,label=label+" fit")
        if label!='': plt.legend()
        err=math.fsum(errfunc(bestfitparams,fitdatax,fitdatay))
        #print 'the best fit has an RMS of {0}'.format(err)
#    plt.t
#    plt.figtext()    
    return bestfitparams
    

def fitlinear(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fits decaying sin wave of form: p[0]*np.sin(2.*pi*p[1]*x+p[2]*pi/180.)*np.e**(-1.*(x-p[5])/p[3])+p[4]"""
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fitting paramters
    if fitparams is None:
        fitparams    = [1,1]
        fitparams[0] = fitdatay[0]
        fitparams[1] = (float(fitdatay[-1])-float(fitdatay[0]))/( float(fitdatax[-1])-float(fitdatax[0]))


    p1 = fitgeneral(fitdatax, fitdatay, linear, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    return p1


def fitpoly(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label="",debug=False):
    """fit lorentzian:
        returns [offset,amplitude,center,hwhm]"""
        
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fitting paramters
    if fitparams is None:
        fitparams=[ydata[0],(ydata[-1]-ydata[0])/(xdata[-1]-xdata[0]),0,xdata[0]]
    
    # fit debugging option : inspection mode
    if debug==True: print(fitparams)
    
    p1 = fitgeneral(fitdatax, fitdatay, poly, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)

    return p1



def fitlor(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label="",debug=False):
    """fit lorentzian:
        returns [offset,amplitude,center,hwhm]"""
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fitting parameters
    if fitparams is None:
        fitparams=[0,0,0,0]
        fitparams[0]=(fitdatay[0]+fitdatay[-1])/2.
        fitparams[1]=max(fitdatay)-min(fitdatay)
        fitparams[2]=fitdatax[np.argmax(fitdatay)]
        fitparams[3]=(max(fitdatax)-min(fitdatax))/10.
    
    # fit debugging option : inspection mode
    if debug==True: print(fitparams)
    
    p1 = fitgeneral(fitdatax, fitdatay, lorfunc, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    p1[3]=abs(p1[3])
    
    return p1



def fitharm(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label="",debug=False):
    """fit lorentzian:
        returns [offset,amplitude,center,hwhm]"""
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fitting parameters
    if fitparams is None:
        fitparams    = [0,0,0,0]
        fitparams[0] = (fitdatay[0]+fitdatay[-1])/2.
        fitparams[3] = (max(fitdatax)-min(fitdatax))/50.
        fitparams[2] = fitdatax[np.argmax(fitdatay)]
        fitparams[1] = (max(fitdatay)-min(fitdatay))*fitparams[3]*fitparams[2]*4*(3.14)**2

#         fitparams[3]=(max(fitdatax)-min(fitdatax))/50.
#         fitparams[3] = 2.88606749e+05
    
    # fit debugging option : inspection mode
    if debug==True: print(fitparams)
    
    p1 = fitgeneral(fitdatax,fitdatay,harmfunc,fitparams,domain=None,showfit=showfit,showstartfit=showstartfit,label=label)
    p1[3] = abs(p1[3])
    p1[2] = abs(p1[2])
    
    return p1



    
def fitexp(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fit exponential decay (p[0]+p[1]*exp(-(x-p[2])/p[3]))"""
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fit parameters
    if fitparams is None:    
        fitparams=[0.,0.,0.,0.]
        fitparams[0]=fitdatay[-1]
        fitparams[1]=fitdatay[0]-fitdatay[-1]
        fitparams[1]=fitdatay[0]-fitdatay[-1]
        fitparams[2]=fitdatax[0]
        fitparams[3]=(fitdatax[-1]-fitdatax[0])/5.
    
    #print fitparams
    p1 = fitgeneral(fitdatax, fitdatay, expfunc, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    return p1   

# Test fit code




def fitpulse_err(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fit pulse err decay (p[0]+p[1]*(1-p[2])^x)"""
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
        
    # determine initial fit parameters
    if fitparams is None:    
        fitparams=[0.,0.]
        fitparams[0]=fitdatay[-1]
        fitparams[1]=fitdatay[0]-fitdatay[-1]
        fitparams[1]=fitdatay[0]-fitdatay[-1]
        
    #print fitparams
    p1 = fitgeneral(fitdatax, fitdatay, pulse_errfunc, fitparams, domain=None, showfit=showfit,
                    showstartfit=showstartfit, label=label)
    return p1   
    



def fitgauss(xdata,ydata,fitparams=None,no_offset=False,domain=None,showfit=False,showstartfit=False,label=""):
    """
    no_offset = True:   p[1] exp(- (x-p[2])**2/p[3]**2/2)
    no_offset = False:  p[0]+p[1] exp(- (x-p[2])**2/p[3]**2/2)
    """
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fit parameters
    if fitparams is None:
        fitparams=[0,0,0,0]
        fitparams[0]=(fitdatay[0]+fitdatay[-1])/2.
        fitparams[1]=max(fitdatay)-min(fitdatay)
        fitparams[2]=fitdatax[np.argmax(fitdatay)]
        fitparams[3]=(max(fitdatax)-min(fitdatax))/3.

    # determine fitfunc (with offset / wo offset)
    if no_offset:
        fitfunc = gaussfunc_nooffset
        fitparams = fitparams[1:]
    else:
        fitfunc = gaussfunc

    p1 = fitgeneral(fitdatax,fitdatay,fitfunc,fitparams,domain=None,showfit=showfit,showstartfit=showstartfit,label=label)
    
    return p1   
    


def fitdecaysin(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fits decaying sin wave of form: p[0]*np.sin(2.*pi*p[1]*x+p[2]*pi/180.)*np.e**(-1.*(x-p[5])/p[3])+p[4]"""
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fit parameters
    if fitparams is None:    
        FFT=scipy.fft.fft(fitdatay)
        fft_freqs=scipy.fftpack.fftfreq(len(fitdatay),fitdatax[1]-fitdatax[0])
        max_ind=np.argmax(abs(FFT[4:int(len(fitdatay)/2)]))+4
        fft_val=FFT[max_ind]
        
        fitparams=[0,0,0,0,0]
        fitparams[4]=np.mean(fitdatay)
        fitparams[0]=(max(fitdatay)-min(fitdatay))/2.#2*abs(fft_val)/len(fitdatay)
        fitparams[1]=fft_freqs[max_ind]
        fitparams[2]=(cmath.phase(fft_val)-np.pi/2.)*180./np.pi
        fitparams[3]=(max(fitdatax)-min(fitdatax))

        #fitparams[5]=fitdatax[0]
        
    decaysin3 = lambda p, x: p[0] * np.sin(2. * np.pi * p[1] * x + p[2] * np.pi / 180.) * np.e ** (
    -1. * (x - fitdatax[0]) / p[3]) + p[4]
    # decaysin3 = lambda p, x: p[0] * np.sin(2. * np.pi * p[1] * x + p[2] - np.pi / 2.) * np.e ** (
    # -1. * (x - fitdatax[0]) / p[3]) + p[4]
    #print "fitparams: ",fitparams
    p1 = fitgeneral(fitdatax, fitdatay, decaysin3, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    return p1  

def fitdecaydoublesin(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fits decaying sin wave of form: p[0]*np.sin(2.*pi*p[1]*x+p[2]*pi/180.)*np.e**(-1.*(x-p[5])/p[3])+p[4]"""
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fit parameters
    if fitparams is None:
        FFT=scipy.fft.fft(fitdatay)
        fft_freqs=scipy.fftpack.fftfreq(len(fitdatay),fitdatax[1]-fitdatax[0])
        max_ind=np.argmax(abs(FFT[4:len(fitdatay)/2.]))+4
        fft_val=FFT[max_ind]

        fitparams=[0,0,0,0,0,0,0,0]
        fitparams[4]=np.mean(fitdatay)
        fitparams[0]=(max(fitdatay)-min(fitdatay))/2.#2*abs(fft_val)/len(fitdatay)
        fitparams[1]=fft_freqs[max_ind]
        fitparams[6]=fft_freqs[max_ind]-0.001
        fitparams[2]=(cmath.phase(fft_val)-np.pi/2.)*180./np.pi
        fitparams[3]=(max(fitdatax)-min(fitdatax))
        fitparams[5] = fitparams[0]
        #fitparams[5]=fitdatax[0]

    decaydoublesin3 = lambda p, x: p[0] * (np.sin(2. * np.pi * p[1] * x + p[2] * np.pi / 180.) + p[5]* np.sin(2. * np.pi * p[6] * x + p[7] * np.pi / 180.) )* np.e ** (
    -1. * (x - fitdatax[0]) / p[3]) + p[4]
    #print "fitparams: ",fitparams
    p1 = fitgeneral(fitdatax, fitdatay, decaydoublesin3, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    return p1


def fitsin(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fits sin wave of form: p[0]*np.sin(2.*pi*p[1]*x+p[2]*pi/180.)+p[3]"""
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fit parameters  
    if fitparams is None:    
        FFT=scipy.fft.fft(fitdatay)
        fft_freqs=scipy.fftpack.fftfreq(len(fitdatay),fitdatax[1]-fitdatax[0])
        max_ind=np.argmax(abs(FFT[4:len(fitdatay)/2.]))+4
        fft_val=FFT[max_ind]
        
        fitparams=[0,0,0,0]
        fitparams[3]=np.mean(fitdatay)
        fitparams[0]=(max(fitdatay)-min(fitdatay))/2.#2*abs(fft_val)/len(fitdatay)
        fitparams[1]=fft_freqs[max_ind]
        fitparams[2]=(cmath.phase(fft_val)-np.pi/2.)*180./np.pi
        #fitparams[3]=(max(fitdatax)-min(fitdatax))
        #fitparams[5]=fitdatax[0]
        
    sin2=lambda p,x: p[0]*np.sin(2.*np.pi*p[1]*x+p[2]*np.pi/180.)+p[3]
    #print "fitparams: ",fitparams
    p1 = fitgeneral(fitdatax, fitdatay, sin2, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    return p1  



    
def fithanger_old (xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fit's Hanger Transmission (S21) data without taking into account asymmetry
       returns p=[f0,Q,S21Min,Tmax]
    """
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fit parameters
    if fitparams is None:    
        fitparams    = [0,0,0,0]
        peakloc      = np.argmin(fitdatay)                  # locate valley 
        ymax         = (fitdatay[0]+fitdatay[-1])/2.        # y max = avg of yvalue at each end
        fitparams[0] = fitdatax[peakloc]                    # set valley as a resonance
        fitparams[1] = abs(fitdatax[peakloc]/((max(fitdatax)-min(fitdatax))/3.))   # Q ~ f0/FWHM
        
        if fitdatay[peakloc] > 0:
            fitparams[2] = (fitdatay[peakloc] / ymax) ** 0.5
        else:
            fitparams[2] = 0.001
        
        fitparams[3] = ymax
   
    return fitgeneral(fitdatax, fitdatay, hangerfunc_old, fitparams, domain=None, showfit=showfit,
                      showstartfit=showstartfit, label=label)




def fithanger_new(xdata, ydata, fitparams=None, domain=None, showfit=False, showstartfit=False, printresult=False,
                  label="", mark_data='bo', mark_fit='r-'):
    """Fit Hanger Transmission (S21) data taking into account asymmetry.
        needs a given Qc, which is assumed to be constant.
        You need to define the Qc in hangerfunc_new()
        fitparams = []
        returns p=[f0,Qi,df,scale]
        Uses hangerfunc_new. 
    """
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax = xdata
        fitdatay = ydata
    
    # determine initial fit parameters
    if fitparams is None:    
        peakloc = np.argmin(fitdatay)
        ymax    = (fitdatay[0]+fitdatay[-1])/2.
        ymin    = fitdatay[peakloc]        
        f0      = fitdatax[peakloc]
        Q0      = abs(fitdatax[peakloc]/((max(fitdatax)-min(fitdatax))/5.))
        scale   = ymax-ymin
        Qi      = 2*Q0

        #slope = (fitdatay[-1]-fitdatay[0])/(fitdatax[-1]-fitdatax[0])
        #offset= ymin-slope*f0
        fitparams = [f0,abs(Qi),0.,scale]
        #print '--------------Initial Parameter Set--------------\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}\nScale: {4}\nSlope: {5}\nOffset:{6}\n'.format(f0,Qi,Qc,0.,scale,slope, offset)
    fitresult    = fitgeneral(fitdatax, fitdatay, hangerfunc_new, fitparams, domain=None, showfit=showfit,
                           showstartfit=showstartfit, label=label, mark_data=mark_data, mark_fit=mark_fit)
    
    fitresult[1] = abs(fitresult[1])
    #fitresult[2]=abs(fitresult[2])
    if printresult: print('-- Fit Result --\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}'.format(fitresult[0], fitresult[1],
                                                                                        fitresult[2], fitresult[3]))
    return fitresult
    

def fithanger_new_withQc(xdata, ydata, fitparams=None, domain=None, showfit=False, showstartfit=False,
                         printresult=False, label="", mark_data='bo', mark_fit='r-'):
    """Fit Hanger Transmission (S21) data taking into account asymmetry.
        use the same parameters as old one 'fithanger', but a different interpretation of the fit formula
        fitparams = []
        returns p=[f0,Qi,Qc,df,scale]
        Uses hangerfunc. 
    """
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
        
    # determine initial fit parameters
    if fitparams is None:    
        peakloc = np.argmin(fitdatay)
        ymax    = (fitdatay[0]+fitdatay[-1])/2.
        ymin    = fitdatay[peakloc]        
        f0      = fitdatax[peakloc]
        Q0      = abs(fitdatax[peakloc]/((max(fitdatax)-min(fitdatax))/5.))
        scale   = ymax-ymin
        Qi      = 2*Q0
        Qc      = Q0
        #slope = (fitdatay[-1]-fitdatay[0])/(fitdatax[-1]-fitdatax[0])
        #offset= ymin-slope*f0
        fitparams = [f0,abs(Qi),abs(Qc),0.,scale]
        #print '--------------Initial Parameter Set--------------\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}\nScale: {4}\nSlope: {5}\nOffset:{6}\n'.format(f0,Qi,Qc,0.,scale,slope, offset)
    
    fitresult = fitgeneral(fitdatax, fitdatay, hangerfunc_new_withQc, fitparams, domain=None, showfit=showfit,
                           showstartfit=showstartfit, label=label, mark_data=mark_data, mark_fit=mark_fit)
    
    fitresult[1] = abs(fitresult[1])
    fitresult[2] = abs(fitresult[2])
    
    if printresult: print('-- Fit Result --\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}\nscale: {4}'.format(fitresult[0],
                                                                                                    fitresult[1],
                                                                                                    fitresult[2],
                                                                                                    fitresult[3],
                                                                                                    fitresult[4]))
    return fitresult

def fitasym_lorentzian(xdata, ydata, fitparams=None, domain=None, showfit=False, showstartfit=False,
                         printresult=False, label="", mark_data='.', mark_fit='r-'):
    
    '''
    data_folder : working directory
    file: filename
    guess : [freq, df, external Q, total Q, scale1, scale2]
    span1 : fit region left
    span2 : fit region right
    
    '''
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
        
    # determine initial fit parameters
    if fitparams is None:    
        peakloc = np.argmin(fitdatay)
        ymax    = (fitdatay[0]+fitdatay[-1])/2.
        ymin    = fitdatay[peakloc]        
        f0      = fitdatax[peakloc]
        Q0      = abs(fitdatax[peakloc]/((max(fitdatax)-min(fitdatax))/5.))
        scale   = ymax-ymin
        Q       = 3/4*Q0
        Qc      = Q0
        #slope = (fitdatay[-1]-fitdatay[0])/(fitdatax[-1]-fitdatax[0])
        #offset= ymin-slope*f0
        fitparams = [f0,0,abs(Qc),abs(Qc),scale]
        
    fitresult = fitgeneral(fitdatax, 10**(fitdatay/20), lorentzian_asym, fitparams, domain=None, showfit=showfit,
                           showstartfit=showstartfit, label=label, mark_data=mark_data, mark_fit=mark_fit)
    
    
    #fitresult[1] = abs(fitresult[1])
    #fitresult[2] = abs(fitresult[2])
    Q_loss = (fitresult[2]*fitresult[3])/(fitresult[2]-fitresult[3])
    if showfit:
        plt.plot(fitresult[0], lorentzian_asym(fitresult, fitresult[0]), 'go', markersize=15, label='fitted coupled Q = {}, \n fitted total Q ={} \n reson freq = {}GHz'.format(fitresult[2],
                                                                                                                                                                                fitresult[3],
                                                                                                                                                                                fitresult[0]))
        plt.xlabel("frequency(GHz)")
        plt.ylabel("S21(lin mag)")
        plt.legend(loc='lower left')
    if printresult: 
        print(fitresult)
        print('-- Fit Result --\nf0: {0}\ndf: {1}\nQc: {2}\nQloss: {3}\nQtot: {4}\nscale: {5}'.format(fitresult[0],
                                                                                                    fitresult[1],
                                                                                                    fitresult[2],
                                                                                                    Q_loss,
                                                                                                    fitresult[3],
                                                                                                    fitresult[4]))
    return fitresult
        

def fithanger(xdata, ydata, fitparams=None, domain=None, showfit=False, showstartfit=False, printresult=False, label="",
              mark_data='bo', mark_fit='r-'):
    """Fit Hanger Transmission (S21) data taking into account asymmetry.
        fitparams = []
        returns p=[f0,Qi,Qc,df,scale]
        Uses hangerfunc. 
    """
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
        
    if fitparams is None:    
        peakloc=np.argmin(fitdatay)
        ymax=(fitdatay[0]+fitdatay[-1])/2.
        ymin=fitdatay[peakloc]        
        f0=fitdatax[peakloc]
        Q0=abs(fitdatax[peakloc]/((max(fitdatax)-min(fitdatax))/3.))
        scale= ymax
        Qi=Q0*(1.+ymax)
        Qc=Qi/(ymax)
        #slope = (fitdatay[-1]-fitdatay[0])/(fitdatax[-1]-fitdatax[0])
        #offset= ymin-slope*f0
        fitparams=[f0,abs(Qi),abs(Qc),0.,scale]
        #print '--------------Initial Parameter Set--------------\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}\nScale: {4}\nSlope: {5}\nOffset:{6}\n'.format(f0,Qi,Qc,0.,scale,slope, offset)
    fitresult = fitgeneral(fitdatax, fitdatay, hangerfunc, fitparams, domain=None, showfit=showfit,
                           showstartfit=showstartfit, label=label, mark_data=mark_data, mark_fit=mark_fit)
    if printresult: print('-- Fit Result --\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}\nScale: {4}'.format(fitresult[0],
                                                                                                    fitresult[1],
                                                                                                    fitresult[2],
                                                                                                    fitresult[3],
                                                                                                    fitresult[4]))
    return fitresult
    

def fithangertilt(xdata, ydata, fitparams=None, domain=None, showfit=False, showstartfit=False, printresult=False,
                  label="", mark_data='bo', mark_fit='r-'):
    """Fit tilted Hanger Transmission (S21) data taking into account asymmetry.
        fitparams = []
        returns p=[f0, Qi, Qc, df, slope, offset]
        Uses hangerfunctilt instead of hangerfunc. 
    """
    
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    
    # determine initial fit parameters
    if fitparams is None:    
        peakloc = np.argmin(fitdatay)
        ymax    = (fitdatay[0]+fitdatay[-1])/2.
        ymin    = fitdatay[peakloc]        
        f0      = fitdatax[peakloc]
        Q0      = abs(fitdatax[peakloc]/((max(fitdatax)-min(fitdatax))/3.))
        Qi      = Q0*(1.+ymax)
        Qc      = Qi/ymax
        scale   = ymax-ymin
        slope   = (fitdatay[-1]-fitdatay[0])/(fitdatax[-1]-fitdatax[0])
        offset  = ymin-slope*f0
        fitparams = [f0,Qi,Qc,0.0001,slope, offset]
        #print '--------------Initial Parameter Set--------------\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}\nScale: {4}\nSlope: {5}\nOffset:{6}\n'.format(f0,Qi,Qc,0.,scale,slope, offset)
        
        fitresult = fitgeneral(fitdatax, fitdatay, hangerfunctilt, fitparams, domain=None, showfit=showfit,
                               showstartfit=showstartfit, label=label)
        
        if printresult: print('-- Fit Result --\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}\nslope: {4}\noffset: {5}\n'.format(
            fitresult[0], fitresult[1], fitresult[2], fitresult[3], fitresult[4], fitresult[5]))
    return fitresult



def fit_SNT(xdata, ydata, fitparams=None, domain=None, showfit=False, showstartfit=False, label='', debug=False):
    """fit Shot Noise Thermometer curve:
        returns [Tn,GainBW,T,voffset]"""
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata

        # Get starting parameters
    if fitparams is None:
        # fit high bias region with linear fit
        # get high bias data
        edge_index = len(xdata) / 3
        edge_x = xdata[-edge_index:]
        A = np.array([edge_x, np.ones(len(edge_x))])
        w = np.linalg.lstsq(A.T, ydata[-edge_index:])[0]

        qe, kb = (1.6e-19, 1.38e-23)
        GB_guess = w[0] / (.5 * qe)
        Tn_guess = w[1] / (kb * GB_guess)
        T_guess = abs((min(ydata) - w[1]) / (kb * GB_guess))
        voff_guess = 0.002
        fitparams = (Tn_guess, GB_guess, T_guess, voff_guess)

    if debug == True: print(fitparams)
    p1 = fitgeneral(fitdatax, fitdatay, SNT_func, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    return p1



def fitbackground(xdata,ydata,fitparams=None, showfit=False,showstartfit=False,label=""):
    """Fit Hanger Transmission (S21) data taking into account asymmetry.
        fitparams = []
        returns p=[f0,Qi,Qc,df,scale]
        Uses hangerfunc. 
    """
    fitdatax=xdata
    fitdatay=ydata
    if fitparams is None:    
        fitparams=[-6,0,0,0,0,0,0,0,0,0,6.9e+9]
        #print '--------------Initial Parameter Set--------------\nf0: {0}\nQi: {1}\nQc: {2}\ndf: {3}\nScale: {4}\nSlope: {5}\nOffset:{6}\n'.format(f0,Qi,Qc,0.,scale,slope, offset)
    return fitgeneral(fitdatax, fitdatay, polynomial, fitparams, domain=None, showfit=showfit,
                      showstartfit=showstartfit, label=label)
     



def fitdecayrabi(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fits decaying sin wave of form: p[0]*np.sin(2.*pi*p[1]*x+p[2]*pi/180.)*np.e**(-1.*(x-p[5])/p[3])+p[4]"""
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    if fitparams is None:
        FFT=scipy.fft.fft(fitdatay)
        fft_freqs=scipy.fftpack.fftfreq(len(fitdatay),fitdatax[1]-fitdatax[0])
        max_ind=np.argmax(abs(FFT[2:int(len(fitdatay)/2.)]))+2
        fft_val=FFT[max_ind]

        fitparams=[0,0,0,0,0,0]
        fitparams[0]= min(fitdatay)
        fitparams[1]= max(fitdatay) - min(fitdatay)
        # print fitparams[1]
        fitparams[2]= 0
        fitparams[3]= (max(fitdatax)-min(fitdatax))
        #fitparams[4]= 2*abs(fft_val)/len(fitdatay)
        fitparams[4]= fft_freqs[max_ind]
        # fitparams[5]= 0  #2*abs(fft_val)/len(fitdatay)

    p1 = fitgeneral(fitdatax,fitdatay,decayrabifunc1,fitparams,domain=None,showfit=showfit,showstartfit=showstartfit,label=label)
    return p1


def fitdoubleexp(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fit exponential decay (p[0]+p[1]*exp(-(x-p[2])/p[3]))"""
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    if fitparams is None:
        fitparams=[0.,0.,0.,0.,0.,0.]
       # fitparams[0]=fitdatay[-1]
        fitparams[1]= max(fitdatay) - fitdatay[0]
        #fitparams[1]=fitdatay[0]-fitdatay[-1]
        #fitparams[2]=fitdatax[0]
        fitparams[3]=(fitdatax[-1]-fitdatax[0])/5
        #fitparams[4]=(max(fitdatay)-fitdatay[-1])
        fitparams[5]=(fitdatax[-1]-fitdatax[0])/5
    #print fitparams
    p1 = fitgeneral(fitdatax,fitdatay,doubleexpfunc,fitparams,domain=None,showfit=showfit,showstartfit=showstartfit,label=label)
    return p1


def fitsin(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label=""):
    """Fits decaying sin wave of form: p[0]*np.sin(2.*pi*p[1]*x+p[2]*pi/180.)*np.e**(-1.*(x-p[5])/p[3])+p[4]"""
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    if fitparams is None:
        FFT=scipy.fft.fft(fitdatay)
        fft_freqs=scipy.fftpack.fftfreq(len(fitdatay),fitdatax[1]-fitdatax[0])
        max_ind=np.argmax(abs(FFT[4:int(len(fitdatay)/2.)]))+4
        fft_val=FFT[max_ind]

        fitparams=[0,0,0,0]
        fitparams[3]=np.mean(fitdatay)
        fitparams[0]=(max(fitdatay)-min(fitdatay))/2.#2*abs(fft_val)/len(fitdatay)
        fitparams[1]=fft_freqs[max_ind]
        fitparams[2]=(cmath.phase(fft_val)-np.pi/2.)*180./np.pi
        fitparams[3]=(max(fitdatax)-min(fitdatax))

        #fitparams[5]=fitdatax[0]

    sin3 = lambda p, x: p[0] * np.sin(2. * np.pi * p[1] * x + p[2] * np.pi / 180.) + p[3]
    #print "fitparams: ",fitparams
    p1 = fitgeneral(fitdatax, fitdatay, sin3, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    return p1





def fitrabisatfunc(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label="",debug=False):
    """fit lorentzian:
        returns [offset,amplitude,center,hwhm]"""
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    if fitparams is None:
        fitparams=[0,1]
    if debug==True: print(fitparams)
    p1 = fitgeneral(fitdatax, fitdatay, rabisatfunc, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)

    return p1


def fitrabiwidth(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label="",debug=False):
    """fit lorentzian:
        returns [offset,amplitude,center,hwhm]"""
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    if fitparams is None:
        fitparams=[np.sqrt((ydata[-1]-ydata[0])/(xdata[-1]-xdata[0])),ydata[0]]
    if debug==True: print(fitparams)
    p1 = fitgeneral(fitdatax, fitdatay, rabiwidth, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)

    return p1

def fitrabiwidth2(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label="",debug=False):
    """fit lorentzian:
        returns [offset,amplitude,center,hwhm]"""
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    if fitparams is None:
        fitparams=[np.sqrt((ydata[-1]-ydata[0])/(xdata[-1]-xdata[0])),ydata[0]]
    if debug==True: print(fitparams)
    p1 = fitgeneral(fitdatax, fitdatay, rabiwidth2, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)

    return p1




def fitdispersiveshift(xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label="",debug=False):
    """fit lorentzian:
        returns [offset,amplitude,center,hwhm]"""
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata
    if fitparams is None:
        fitparams=[0,0,0,0,0,0,0]
        fitparams[0]=(fitdatay[0]+fitdatay[-1])/2.
        fitparams[1]=max(fitdatay)-min(fitdatay)
        fitparams[2]=fitdatax[np.argmax(fitdatay)]
        fitparams[3]=(max(fitdatax)-min(fitdatax))/10.
        fitparams[4]=0
        fitparams[5]=fitdatax[np.argmax(fitdatay)]
        fitparams[6]=(max(fitdatax)-min(fitdatax))/10.
    if debug==True: print(fitparams)
    p1 = fitgeneral(fitdatax, fitdatay, dispersiveshift, fitparams, domain=None, showfit=showfit, showstartfit=showstartfit,
                    label=label)
    p1[3]=abs(p1[3])
    return p1


def fitshortABSspectrum(xdata, ydata, fitparams=None, domain=None, showfit=False, showstartfit=False, label="",debug=False):
    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata
        fitdatay=ydata

def fitTLSmodel(f, T, tls_sat_point, xdata,ydata,fitparams=None,domain=None,showfit=False,showstartfit=False,label="",debug=False):

    def specific_tlsmodel(p, x):
        p0, p1, p2, p3 = p
        y = tlsfit_model([p0, p1, p2, p3, f, T], x)
        return y

    if domain is not None:
        fitdatax,fitdatay = selectdomain(xdata,ydata,domain)
    else:
        fitdatax=xdata[:tls_sat_point]   # e.g. power_circ[:-3]
        fitdatay=ydata[:tls_sat_point]   # e.g.  1/Qint[:-3]
    if fitparams is None:
        print("Assign initial guess... this fitting is unstable")
        
    if debug==True: print(fitparams)

    TLSfit = fitgeneral(xdata=fitdatax, ydata=fitdatay, fitfunc=specific_tlsmodel, fitparams=[1e-4, 1e-5, 0.01, 1e-12], domain=None, showfit=True, showstartfit=False, showdata=True, label='', mark_data='bo', mark_fit='r-')
    print(f"Q0 = {1/TLSfit[0]}, losstan = {TLSfit[1]}, beta = {TLSfit[2]}, Pc = {TLSfit[3]},")
    plt.ylabel(r"$(Q_{{int}})^{{-1}}$")
    plt.xlabel("circulating power (W)")
    return TLSfit

# In[debug]:
    
if __name__ =='__main__':
    plt.figure(1)
    xdata=np.linspace(-15,25,1000)
    
    params=[1.,20.,5.,2.]
    ydata=gaussfunc(params,xdata)-1+2*np.random.rand(len(xdata))
    #plot(xdata,ydata-2.5+5*random.rand(xdata.__len__()),'bo')
    plt.subplot(1,2,1)
    p1=fitgauss(xdata,ydata,showfit=True)
    plt.subplot(1,2,2)
    p2=fitlor(xdata,ydata,showfit=True)
    #plot(xdata,lorfunc(p1,xdata),'r-')
    
    
    noise=0.
    plt.figure(2)
    params2=[7.8,200,200.,0.005,1.,0.,0.]
    print('{0}\n--------------Test Parameter---------- \nf0: {1}\nQi: {2}\nQc: {3}\ndf: {4}\nScale: {5}\nSlope: {6}\nOffset:{7}\n'.format\
          ('',params2[0],params2[1],params2[2],params2[3],params2[4],params2[5],params2[6]))
#    params2=[7.8,200,0.01,1.]
    xdata2=np.linspace(7,9,1000)
    ydata2=hangerfunc(params2,xdata2)-noise/2.+noise*np.random.rand(len(xdata2))
    fit=fithanger(xdata2,ydata2,showfit=True,showstartfit=True)   
    print('{0}\n--------------Best Fit---------- \nf0: {1}\nQi: {2}\nQc: {3}\ndf: {4}\nScale: {5}\nSlope: {6}\nOffset:{7}\n'.format\
          ('hanger',fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6]))
    #print hangerqs(p3)
    
    plt.show()




# In[6]: peak detect

def _datacheck_peakdetect(x_axis, y_axis):
    if x_axis is None:
        x_axis = list(range(len(y_axis)))
    
    if len(y_axis) != len(x_axis):
        raise ValueError
    
    #needs to be a numpy array
    y_axis = np.array(y_axis)
    x_axis = np.array(x_axis)
    return x_axis, y_axis
    
def peakdetect(y_axis, x_axis = None, lookahead = 300, delta=0):
    """
    Converted from/based on a MATLAB script at: 
    http://billauer.co.il/peakdet.html
    
    function for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively
    
    keyword arguments:
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- (optional) A x-axis whose values correspond to the y_axis list
        and is used in the return to specify the postion of the peaks. If
        omitted an index of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 200) 
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the function from picking up false peaks towards to end of
        the signal. To work well delta should be set to delta >= RMSnoise * 5.
        (default: 0)
            delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the function
    
    return -- two lists [max_peaks, min_peaks] containing the positive and
        negative peaks respectively. Each cell of the lists contains a tupple
        of: (position, peak_value) 
        to get the average peak value do: np.mean(max_peaks, 0)[1] on the
        results to unpack one of the lists into x, y coordinates do: 
        x, y = zip(*tab)
    """
    max_peaks = []
    min_peaks = []
    dump = []   #Used to pop the first hit which almost always is false
       
    # check input data
    x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis)
    # store data length for later use
    length = len(y_axis)
    
    
    #perform some checks
    if lookahead < 1:
        raise ValueError("Lookahead must be '1' or above in value")
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError("delta must be a positive number")
    
    #maxima and minima candidates are temporarily stored in
    #mx and mn respectively
    mn, mx = np.Inf, -np.Inf
    
    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], 
                                        y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x
        
        ####look for max####
        if y < mx-delta and mx != np.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue
            #else:  #slows shit down this does
            #    mx = ahead
            #    mxpos = x_axis[np.where(y_axis[index:index+lookahead]==mx)]
        
        ####look for min####
        if y > mn+delta and mn != -np.Inf:
            #Minima peak candidate found 
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                min_peaks.append([mnpos, mn])
                dump.append(False)
                #set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
            #else:  #slows shit down this does
            #    mn = ahead
            #    mnpos = x_axis[np.where(y_axis[index:index+lookahead]==mn)]
    
    
    #Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            max_peaks.pop(0)
        else:
            min_peaks.pop(0)
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass
        
    return [max_peaks, min_peaks]