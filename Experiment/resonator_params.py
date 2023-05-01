# In[0]: necessary packages

import numpy as np
from numpy import pi, cosh, sinh, tanh, sqrt, log10
from scipy.special import ellipk

# In[1]: functions

# Constants

epsilon_0 = 8.85e-12
k_b=1.380649e-23
h = 6.626068e-34 
c=3e8
mu=4*pi*(1e-7)

'''
take values of ABCD and characteristic impedance and return S-matrix
'''

def abcd2s(A,B,C,D,z0):
    divider = A+(B/z0)+(C*z0)+D
    s11new = (A+(B/z0)-(C*z0)-D)/divider
    s12new = (2*(A*D-B*C))/divider
    s21new = 2/divider
    s22new = (-A+(B/z0)-(C*z0)+D)/divider
    s_new = np.array([[s11new, s12new], [s21new, s22new]])
    return s_new

'''
CPWabcd(L, frequency, impedance, epseff, attenuationConst)
 L = length -> meters
 frequency -> Hz
 impedance -> ohms
 epseff -> unitless (effective dielectric constant)
 attenuationConst -> nepers/m (dB/m/8.686 = Nepers/m)
'''
def coupling_abcd(frequency, Cg, Lg, coupling):
    if coupling == 'L':
        Z = 1j*2*pi*frequency*Lg/2
    elif coupling =='C':
        Z = 2/(1j*2*pi*frequency*Cg)
        
    abcd = np.array([[1,Z],[0,1]])
    return abcd

def reson_abcd(frequency, Cg, f0, R):
    Y = 1j*Cg*2*pi*frequency/(1-(frequency/f0)**2+1j*2*pi*frequency*R*Cg)
    abcd = np.array([[1,0],[Y,1]])
    return abcd
        
def CPWabcd(L, frequency, impedance, epseff, attenuationConst):

    lambda0 = 3e8/frequency # wavelength
    beta0 = 2*pi/lambda0
    alpha = attenuationConst # attenuation constant in nepers/m to get loss in dB/length = Nepers/length * 8.686
    beta = beta0*sqrt(epseff) # propogation constant
    gamma = alpha + 1j*beta
    
    Zc = impedance  # characteristic impedance of transmission line;
    a = cosh(gamma*L)
    b = Zc*sinh(gamma*L)
    c = 1/Zc * sinh(gamma*L)
    d = cosh(gamma*L)
    abcd = np.array([[a, b], [c, d]])
    return abcd

def epsilon_eff(a,b,height,erel):
    k = a/b
    k3 = tanh(pi*a/4/height)/tanh(pi*b/4/height)
    kef = ellipk(sqrt(1-k**2))*ellipk(k3) / (ellipk(k) * ellipk(sqrt(1-k3**2)))
    eff = (1+erel*kef)/(1+kef)
    Z = (60*pi/sqrt(eff))/(ellipk(k)/ellipk(sqrt(1-k**2))+ellipk(k3)//ellipk(1-k3**2))
    return (eff,Z)

def L_kinetic_v2(r, tc, width, thick):
    lambda0 = 1.05e-3*sqrt(r/tc)
    Lk = ((mu*lambda0*log10((4*width)/thick)*sinh(thick/lambda0))/(width*(pi**2)*(cosh(thick/lambda0)-1)))
    return (Lk, Lk*width)

def L_kinetic_v3(rho, tc, width, thick, temp):
    D = 1.76*k_b*tc
    Rsq = rho/thick
    Lk=((Rsq*h)/(2*(pi**2)*D*tanh(D/(2*k_b*temp))))/(width)  # per unit length
    return (Lk, Lk*width)

def impe_char(Erel, gap, rho, tc, width, thick, temp):
    epseff=(Erel+1)/2
    k1=width/(width+(2*gap))
    k2=sqrt(1-(k1**2))
    Lg = (mu*ellipk(k2))/(4*ellipk(k1))
    C=(4*epsilon_0*epseff*ellipk(k1))/(ellipk(k2))
    Lk = L_kinetic_v3(rho, tc, width, thick, temp)[0]
    Z=sqrt((Lg+Lk)/C)
    #f_res=1/((sqrt((L_g+L_k3)*C))*(2*l_r))
    
    return Z

#L_kinetic_v3(r=1.5e-7, tc=8, width=10e-6, thick=100e-9, temp=1)

'''
Jiansong Gao Thesis
'''

def k(a, b, thick):
    d=2*thick/pi
    if thick > 0:
        u1=a/2+d/2+3*d/2*log10(2)-d/2*log10(2*d/a)+d/2*log10((b-a)/(b+a))
        u1p=u1-d
        u2=b/2-d/2-3*d/2*log10(2)+d/2*log10(2*d/b)+d/2*log10((b-a)/(b+a))
        u2p=u2+d
    elif thick == 0:
        u1 = a/2
        u2 = b/2
    kt = u1/u2
    ktp = sqrt(1-kt**2)
    return (kt, ktp)

def g(a, b, thick):
    k = a/b
    g0 = 1/(4*a*ellipk(k)**2*(1-k**2))
    g_ctr = g0 * (pi+log10(4*pi*a/thick)-k*log10((1+k)/(1-k)))
    g_gnd = g0 * k *(pi+log10(4*pi*b/thick)-k*log10((1+k)/(1-k)))
    gtot = g_ctr+g_gnd
    return (gtot, g_ctr, g_gnd)

def Cap(a, b, thick):
    eps0 = 8.85e-12
    kt = k(a,b,thick)[0]
    ktp = k(a,b,thick)[1]
    C_half = eps0*2*ellipk(kt)/ellipk(ktp)
    return C_half

def Ind(a, b, thick):
    mu0 = 4*pi*(1e-7)
    kt = k(a,b,thick)[0]
    ktp = k(a,b,thick)[1]
    L_half = mu0*ellipk(ktp)/ellipk(kt)/2
    return L_half

def Ctot(a, b, thick, Erel):
    Ct = Cap(a,b,thick)+Erel*Cap(a,b,0)
    return Ct

def Ceff(a, b, thick, Eeff):
    Ct = Eeff*Cap(a,b,0)
    return Ct
    

def Ltot(a, b, thick):
    Lt = 0.5*Ind(a,b,thick)
    return Lt

def imp(a,b, thick, rho, tc, Erel, temp):
    Lg = Ltot(a,b, thick)
    Lk = L_kinetic_v3(rho, tc, a, thick, temp)[0]
    C = Ctot(a,b, thick, Erel)
    Z = sqrt((Lg+Lk)/C)
    v_ph = 1/sqrt((Lg+Lk)*C)
    return Z, v_ph

def vph(l ,L ,C, QWR):
    num = len(l)
    sumL = 0
    sumC = 0
    for i in range(num):
        sumL += l[i]*L[i]
        sumC += l[i]*C[i]
    totlen = np.sum(l)
    l = sumL/totlen ; c = sumC/totlen
    v= 1/sqrt(l*c)
    wavelen = totlen*4/QWR
    f = v/wavelen
    return (v, wavelen, f) 