# In[0]:
import numpy as np
from numpy import cos, sin, tan, arctan, arcsin, arccos, pi, cosh, sinh, tanh, sqrt, log10
from scipy.special import ellipk
import matplotlib.pyplot as plt
from pyautocad import Autocad, APoint, aDouble
import sys
sys.path.append('C:/Users/20141/OneDrive - purdue.edu/PhD_research/resonators/python_codes/Experiment/')
import os

# In[1]: cpw design parameters 

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

# In[2]: cpw design functions

acad = Autocad(create_if_not_exists=True)

def cpw_linear(a, b, l, init_p, axis, angle=0):
    width = 0.5*(b-a)
    
    
    if axis == 0:  # x-axis
    
        st = [init_p[0]-a/2*sin(angle), init_p[1]+a/2*cos(angle)]
    
        st.append(st[-2]+l*cos(angle))
        st.append(st[-2]+l*sin(angle))
        st.append(st[-2]-width*sin(angle))
        st.append(st[-2]+width*cos(angle))
        st.append(st[-2]-l*cos(angle))
        st.append(st[-2]-l*sin(angle))
        st.append(st[0])
        st.append(st[1])
        
        sq = aDouble(st)
    
        acad.model.AddLightWeightPolyline(sq)
            
        for i in range(int(len(sq)/2)):
            sq[2*i] = sq[2*i]+(a+b)/2*sin(angle)
            sq[2*i+1] = sq[2*i+1]-(a+b)/2*cos(angle)
        
        acad.model.AddLightWeightPolyline(sq) 
    
        terminate = [init_p[0]+l*cos(angle), init_p[1]+l*sin(angle)]
        
    elif axis == 1:  # y-axis
    
        st = [init_p[0]+a/2, init_p[1]]
        
        st.append(st[-2])
        st.append(st[-2]+l)
        st.append(st[-2]+width)
        st.append(st[-2])
        st.append(st[-2])
        st.append(st[-2]-l)
        st.append(st[0])
        st.append(st[1])
        
        sq = aDouble(st)
        
        acad.model.AddLightWeightPolyline(sq)
            
        for i in range(int(len(sq)/2)):
            sq[2*i] = sq[2*i]-(a+b)/2
        
        acad.model.AddLightWeightPolyline(sq) 
    
        terminate = [init_p[0], init_p[1]+l]
        
    return terminate

def polar_sq(origin, azimuth, radius, width, angle):
  
    st = [origin[0]+radius*cos(azimuth), origin[1]+radius*sin(azimuth)]
    st.append(origin[0] + radius*cos(azimuth-angle))
    st.append(origin[1] + radius*sin(azimuth-angle))
    st.append(origin[0] + (radius+width)*cos(azimuth-angle))
    st.append(origin[1] + (radius+width)*sin(azimuth-angle))
    st.append(origin[0] + (radius+width)*cos(azimuth))
    st.append(origin[1] + (radius+width)*sin(azimuth))
    st.append(st[0])
    st.append(st[1])
    
    return st

def cpw_curve(a, b, l, init_p, curve, NOA, axis, circ, init_angle=0):
    # axis 1 : y  axis 0 : x
    # circ 0 : cw  circ 1 :ccw
    theta = curve/NOA
    r = l/theta
    
    
    width = 0.5*(b-a)
    if init_angle != 0:
        origin = [init_p[0]-r*cos(init_angle), init_p[1]-r*sin(init_angle)]
        azi = init_angle
     
    elif axis == 0: 
        origin = [init_p[0], init_p[1]-1*r]
        azi = pi/2
      
    elif axis == 1:
        origin = [init_p[0]+1*r, init_p[1]]
        azi = pi
        
    
    rin = r-b/2
    rout = r+a/2
    
    for i in range(NOA):
    
        plin = polar_sq(origin, azi - (-1)**circ * i*theta, rin, width, (-1)**circ *theta)
        plout = polar_sq(origin, azi - (-1)**circ * i*theta, rout, width, (-1)**circ *theta)
        
        acad.model.AddLightWeightPolyline(aDouble(plin))
        acad.model.AddLightWeightPolyline(aDouble(plout))
        
        terminate = [origin[0]+r*cos(azi-(-1)**circ*curve) , origin[1]+r*sin(azi-(-1)**circ*curve) ]
        
    
    return terminate

def Ind(a,h,N,d,w):
    return 0.00266*a**0.0603*h**0.4429*N**0.954*d**0.606*w**-0.173


# In[2]: neat example [2023 03/16 MoRe 20nm cavity 6GHz on Si]

'''
QWR part 
'''
A = 43
B = 65
rho_MoRe20nm = 1e-6
MoRethick = 20e-9
Tc_MoRe20nm = 8
Erel_Si = 11.7
cavity_mode = 6e9

z = imp(a = A*1e-6, b = B*1e-6, thick = MoRethick, rho = rho_MoRe20nm, tc = Tc_MoRe20nm, Erel=Erel_Si, temp=1)
lk1=L_kinetic_v3(rho = rho_MoRe20nm, tc=Tc_MoRe20nm, width=A*1e-6, thick=MoRethick, temp=1)
lg = Ltot(a = A*1e-6, b = B*1e-6, thick=MoRethick )
cg =  Ctot(a = A*1e-6, b = B*1e-6, thick=MoRethick , Erel=Erel_Si)
v_ph = z[1]
wavelen = v_ph/cavity_mode
l = 4/4*wavelen
print(lk1, lg, cg, l, 1/4*l, z[0])


# In[CAD]
'''
CAD
'''
qwr = l/4*1e6
num_of_angle = 10


print(qwr)

parity=0

# In[drawing - qwr start]
l_start = 1500
point = cpw_linear(a=A, b=B, l=(-1)**parity*l_start, init_p=[0,0], axis=0, angle = pi)
qwr = qwr-l_start

l_ramp = B*pi/4
point = cpw_curve(a=A, b=B, l=l_ramp, init_p=point, curve=-pi, NOA=num_of_angle, axis=0, circ=parity, init_angle=3*pi/2)
#parity = parity+(-1)**parity
qwr = qwr-l_ramp

l_01 =800 #100
point = cpw_linear(a=A, b=B, l=l_01, init_p=point, axis=0, angle=0)
qwr = qwr-l_01
print(qwr)
# In[qwr-meander]
l_01 = 900
meander_lin = 2
while True:
    if qwr >= meander_lin*l_01:
        point = cpw_curve(a=A, b=B, l=l_ramp, init_p=point, curve=pi, NOA=num_of_angle, axis=0, circ=parity)
        qwr = qwr-l_ramp
        parity = parity+(-1)**parity
        print(qwr)
        
        if qwr >= meander_lin*l_01:
        
            point = cpw_linear(a=A, b=B, l=(-1)**parity*l_01, init_p=point, axis=0)
            qwr = qwr - l_01
            print(qwr)
            
            if qwr >= meander_lin*l_01:
        
                point = cpw_curve(a=A, b=B, l=l_ramp, init_p=point, curve=pi, NOA=num_of_angle, axis=0, circ=parity)
                qwr = qwr - l_ramp
                parity = parity+(-1)**parity
                print(qwr)
                
                if qwr >= meander_lin*l_01:
            
                    point = cpw_linear(a=A, b=B, l=(-1)**parity*l_01, init_p=point, axis=0)
                    qwr = qwr - l_01
                    print(qwr)
    else:
        break
    

print(qwr)

point = cpw_curve(a=A, b=B, l=l_ramp, init_p=point, curve=pi, NOA=num_of_angle, axis=0, circ=parity)
qwr = qwr-l_ramp
parity = parity+(-1)**parity
print(qwr)

l_02 = 1.4*l_01

point = cpw_linear(a=A, b=B, l=(-1)**parity*l_02, init_p=point, axis=0)
qwr = qwr-l_02
print(qwr)

point = cpw_curve(a=A, b=B, l=l_ramp, init_p=point, curve=-pi/2, NOA=8, axis=0, circ=parity)
qwr = qwr-l_ramp
parity = parity+(-1)**parity
print(qwr)  

point = cpw_linear(a=A, b=B, l=(-qwr), init_p=point, axis=0, angle=3*pi/2)
qwr = qwr-qwr
print(qwr)