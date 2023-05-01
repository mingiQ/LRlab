# In[0]: necessary packages

import numpy as np
import matplotlib.pyplot as plt
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

def Ltot(a, b, thick):
    Lt = 0.5*Ind(a,b,thick)
    return Lt

def imp(a,b, thick, rho, tc, Erel, temp):
    Lg = Ltot(a,b, thick)
    Lk = L_kinetic_v3(rho, tc, a, thick, temp)[0]
    C = Ctot(a,b, thick, Erel)
    Z = sqrt((Lg+Lk)/C)
    f = 1/sqrt((Lg+Lk)*C)
    return Z, f

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
    

# In[2]:

'''
% PBG simulator
% This program simulates the microwave transmission spectrum for a photonic
% bandgap filter/resonator using the ABCD matrices for coplanar waveguides .


%% Define the CPW paramters
% a period is defined as a length of CPW with impedance 'Z1' followed by a length of CPW with impedance 'Z2' 
'''

lp = 3000e-6 # length of half of a period in m .7000um
lh = 5000e-6 
#f = np.arange(2e9, 12e9, 1e4) # frequency range to simulate plot in Hz 
f = np.arange(2e9, 8e9, 1e4) # frequency range to simulate plot in Hz 

Z1 = imp(a=100e-6, b=110e-6, thick=100e-9, rho=1.1e-6, tc=11.5, Erel=11, temp=1) #31.176 # impedance of the first half period in Ohms 
Z2 = imp(a=15e-6, b=110e-6, thick=100e-9, rho=1.1e-6, tc=11.5, Erel=11, temp=1) #176.075  # impedance of second half period in Ohms
#Cg = 10e-12  # in F (1pF)
#Lg = 1e-3  # in H
#coup = 'C'

epseff1 = 6 #7.389  # effective dielectric constant of substrate
alpha = 0.000/8.686  # microwave loss in nepers/m;

'''
%% Define the PBG defect - use if you want to make a PBG resonator as opposed to a filter .
'''

ld1 = 700e-6  # length of the defect (1/2 lambda resonator) in m .13mm --> 4.5G
ld2 = 5000e-6 
ld3 = 11000e-6
Zd = imp(a=10e-6, b=16e-6, thick=100e-9, rho=1.1e-6, tc=11.5, Erel=11, temp=1) #50.0  # impedance of of the defect (1/2 lambda resonator) in Ohms .
epseff2 = 6  # effective dielectric constant of substrate
alpha = 0.000/8.686  # loss in nepers/m;
s21left = []
s21left2 = [] 
s21right = []
s21right2 = []
s21 = []

for freq in f:
    
    #abcdc = coupling_abcd(freq, Cg, Lg, coup)
    
    abcd1 = CPWabcd(lp, freq, Z1, epseff1, alpha) # construct the ABCD matrix for the length of CPW with impedance 'Z1'
    abcd2 = CPWabcd(lp, freq, Z2, epseff1, alpha) # construct the ABCD matrix for the length of CPW with impedance 'Z2'
    abcdD1 = CPWabcd(ld1, freq, Zd, epseff2, alpha) # construct the ABCD matrix for the length of CPW with impedance 'Zd'
    
    abcd11 = CPWabcd(lh, freq, Z1, epseff1, alpha)
    abcd22 = CPWabcd(lh, freq, Z2, epseff1, alpha)
    abcdD2 = CPWabcd(ld2, freq, Zd, epseff2, alpha)
    
    abcdD3 = CPWabcd(ld3, freq, Zd, epseff2, alpha)
    
    '''
    % The following line constructs the overall ABCD matrix based on the
    % sample geometry . To change the geometry simply change the matrices
    % being multiplied together . For example: to simulate two periods of
    % the CPW bragg reflector followed by a defect followed by one period
    % of the bragg reflector just multipy
    % abcd1*abcd2*abcd1*abcd2*abcdD*abcd1*abcd2
    '''
    #abcdtot1 = abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcdD1@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1
    #abcdtot2 = abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcdD1@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11
    #abcdtot1 = abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcdD1@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcdD2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2
    abcdtot2 = abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcdD1@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11   
    # The next line of code converts the ABCD matrix to S paramters
    #abcdtot = abcdtot1@abcdtot2
    #abcdtot_right = abcdtot1
    abcdtot_left = abcdtot2
    #abcdtot_right2 = abcdtot1@abcdtot1
    #abcdtot_left2 = abcdtot2@abcdtot2
    
    #sp = abcd2s(abcdtot[0,0], abcdtot[0,1], abcdtot[1,0], abcdtot[1,1], 50)
    #sp_right = abcd2s(abcdtot_right[0,0], abcdtot_right[0,1], abcdtot_right[1,0], abcdtot_right[1,1], 50)
    #sp_right2 = abcd2s(abcdtot_right2[0,0], abcdtot_right2[0,1], abcdtot_right2[1,0], abcdtot_right2[1,1], 50)
    sp_left = abcd2s(abcdtot_left[0,0], abcdtot_left[0,1], abcdtot_left[1,0], abcdtot_left[1,1], 50)
    #sp_left2 = abcd2s(abcdtot_left2[0,0], abcdtot_left2[0,1], abcdtot_left2[1,0], abcdtot_left2[1,1], 50)
    
    # The following code simply picks the s paramter you want to plot and
    # saves it to the variable 's'
    s21left.append(sp_left[1,0])
    #s21left2.append(sp_left2[1,0])
    #s21right.append(sp_right[1,0])
    #s21right2.append(sp_right2[1,0])
    #s21.append(sp[1,0])

#plt.plot(f/1e9, 20*np.log10(np.abs(s21))) 
plt.plot(f/1e9, 20*np.log10(np.abs(s21left))) 
#plt.plot(f/1e9, 20*np.log10(np.abs(s21right)))
#plt.plot(f/1e9, 20*np.log10(np.abs(s21right2))) 
plt.xlabel("Frequency(GHz)")
plt.ylabel("S21(dB)")

# In[maxfreq]
S21 = 20*np.log10(np.abs(s21left))
newS21 = S21[np.where((f>4e9) & (f <6e9))[0]]
newf = f[np.where((f>4e9) & (f <6e9))[0]]

# In[3]: plot data
plt.plot(f/1e9, 20*np.log10(np.abs(s21)), label='Bragg1+Bragg2') 
plt.plot(f/1e9, 20*np.log10(np.abs(s21left)),label='Bragg1') 
plt.plot(f/1e9, 20*np.log10(np.abs(s21right)),label='Bragg2') 
plt.xlabel("Frequency(GHz)")
plt.ylabel("S21(dB)")
plt.legend()

# In[3]: plot data
 
plt.plot(f/1e9, 20*np.log10(np.abs(s21left2)),label='Bragg1') 
plt.plot(f/1e9, 20*np.log10(np.abs(s21left)),label='Bragg1+Bragg1') 
plt.xlabel("Frequency(GHz)")
plt.ylabel("S21(dB)")
plt.legend()

# In[4]: function define 

#lp = 3000e-6 # length of half of a period in m .7000um
#lh = 5000e-6 

def PBGsim(lh, ld1, no_of_unit_cells, f_start, f_stop, res):

    f = np.arange(f_start, f_stop, res) # frequency range to simulate plot in Hz 
    
    Z1 = imp(a=100e-6, b=110e-6, thick=100e-9, rho=1.1e-6, tc=11.5, Erel=11, temp=1) #31.176 # impedance of the first half period in Ohms 
    Z2 = imp(a=15e-6, b=110e-6, thick=100e-9, rho=1.1e-6, tc=11.5, Erel=11, temp=1) #176.075  # impedance of second half period in Ohms

    
    epseff1 = 6 #7.389  # effective dielectric constant of substrate
    alpha = 0.000/8.686  # microwave loss in nepers/m;
    
    '''
    %% Define the PBG defect - use if you want to make a PBG resonator as opposed to a filter .
    '''
    
    Zd = imp(a=10e-6, b=16e-6, thick=100e-9, rho=1.1e-6, tc=11.5, Erel=11, temp=1) #50.0  # impedance of of the defect (1/2 lambda resonator) in Ohms .
    epseff2 = 6  # effective dielectric constant of substrate
    alpha = 0.000/8.686  # loss in nepers/m;
    s21left = []
    
    for freq in f:
               
        abcd11 = CPWabcd(lh, freq, Z1, epseff1, alpha)
        abcd22 = CPWabcd(lh, freq, Z2, epseff1, alpha)
        abcdD1 = CPWabcd(ld1, freq, Zd, epseff2, alpha)
        
        '''
        % The following line constructs the overall ABCD matrix based on the
        % sample geometry . To change the geometry simply change the matrices
        % being multiplied together . For example: to simulate two periods of
        % the CPW bragg reflector followed by a defect followed by one period
        % of the bragg reflector just multipy
        % abcd1*abcd2*abcd1*abcd2*abcdD*abcd1*abcd2
        '''
        
        abcdtot2 = abcdD1
        for i in range(no_of_unit_cells):
            abcdtot2 = abcd11@abcd22@abcdtot2
        for j in range(no_of_unit_cells):
            abcdtot2 = abcdtot2@abcd22@abcd11
        
        #abcdtot2 = abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcdD1@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11   
       
        sp_left = abcd2s(abcdtot2[0,0], abcdtot2[0,1], abcdtot2[1,0], abcdtot2[1,1], 50)
        s21left.append(sp_left[1,0])
    
    
    plt.plot(f/1e9, 20*np.log10(np.abs(s21left))) 
    
    plt.xlabel("Frequency(GHz)")
    plt.ylabel("S21(dB)")
    
# In[4.1]: notes

PBGsim(lh=5000e-6, ld1=700e-6, no_of_unit_cells=5, f_start=2e9, f_stop=8e9, res=1e5)
# In[4.1]: notes

PBGsim(lh=6000e-6, ld1=1000e-6, no_of_unit_cells=4, f_start=2e9, f_stop=8e9, res=1e5)
# In[4.1]: notes

PBGsim(lh=6000e-6, ld1=1000e-6, no_of_unit_cells=5, f_start=2e9, f_stop=8e9, res=1e5)

# In[2]:

'''
% PBG simulator
% This program simulates the microwave transmission spectrum for a photonic
% bandgap filter/resonator using the ABCD matrices for coplanar waveguides .


%% Define the CPW paramters
% a period is defined as a length of CPW with impedance 'Z1' followed by a length of CPW with impedance 'Z2' 
'''

lp = 3000e-6 # length of half of a period in m .7000um
lh = 5000e-6 
f = np.arange(2e9, 12e9, 1e4) # frequency range to simulate plot in Hz 

Z1 = 30 # impedance of the first half period in Ohms
Z2 = 120  # impedance of second half period in Ohms
Cg = 10e-12  # in F (1pF)
Lg = 1e-3  # in H
coup = 'C'

epseff1 = 6 #7.389  # effective dielectric constant of substrate
alpha = 0.000/8.686  # microwave loss in nepers/m;

'''
%% Define the PBG defect - use if you want to make a PBG resonator as opposed to a filter .
'''

ld1 = 8000e-6  # length of the defect (1/2 lambda resonator) in m .13mm --> 4.5G
ld2 = 500e-6 
Zd = 50.0  # impedance of of the defect (1/2 lambda resonator) in Ohms .
epseff2 = 6  # effective dielectric constant of substrate
alpha = 0.000/8.686  # loss in nepers/m;
s21left = []
s21left2 = [] 
s21right = []
s21right2 = []
s21 = []

for freq in f:
    
    abcdc = coupling_abcd(freq, Cg, Lg, coup)
    
    abcd1 = CPWabcd(lp, freq, Z1, epseff1, alpha) # construct the ABCD matrix for the length of CPW with impedance 'Z1'
    abcd2 = CPWabcd(lp, freq, Z2, epseff1, alpha) # construct the ABCD matrix for the length of CPW with impedance 'Z2'
    abcdD1 = CPWabcd(ld1, freq, Zd, epseff2, alpha) # construct the ABCD matrix for the length of CPW with impedance 'Zd'
    
    
    abcd11 = CPWabcd(lh, freq, Z1, epseff1, alpha)
    abcd22 = CPWabcd(lh, freq, Z2, epseff1, alpha)
    abcdD2 = CPWabcd(ld2, freq, Zd, epseff2, alpha)
    reson1 = reson_abcd(freq, 1e-9, 3e9, 0.001)
    reson2 = reson_abcd(freq, 1e-9, 7e9, 0.001)
    
    '''
    % The following line constructs the overall ABCD matrix based on the
    % sample geometry . To change the geometry simply change the matrices
    % being multiplied together . For example: to simulate two periods of
    % the CPW bragg reflector followed by a defect followed by one period
    % of the bragg reflector just multipy
    % abcd1*abcd2*abcd1*abcd2*abcdD*abcd1*abcd2
    '''
    
    #abcdtot1 = abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1@abcd2@abcd1
    abcdtot2 = abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@abcd22@abcd11@reson1@reson2
    # The next line of code converts the ABCD matrix to S paramters
    #abcdtot = abcdtot1@abcdtot2
    #abcdtot_right = abcdtot1
    abcdtot_left = abcdtot2
    #abcdtot_right2 = abcdtot1@abcdtot1
    abcdtot_left2 = abcdtot2@abcdtot2
    
    #sp = abcd2s(abcdtot[0,0], abcdtot[0,1], abcdtot[1,0], abcdtot[1,1], 50)
    #sp_right = abcd2s(abcdtot_right[0,0], abcdtot_right[0,1], abcdtot_right[1,0], abcdtot_right[1,1], 50)
    #sp_right2 = abcd2s(abcdtot_right2[0,0], abcdtot_right2[0,1], abcdtot_right2[1,0], abcdtot_right2[1,1], 50)
    sp_left = abcd2s(abcdtot_left[0,0], abcdtot_left[0,1], abcdtot_left[1,0], abcdtot_left[1,1], 50)
    sp_left2 = abcd2s(abcdtot_left2[0,0], abcdtot_left2[0,1], abcdtot_left2[1,0], abcdtot_left2[1,1], 50)
    
    # The following code simply picks the s paramter you want to plot and
    # saves it to the variable 's'
    s21left.append(sp_left[1,0])
    s21left2.append(sp_left2[1,0])
    #s21right.append(sp_right[1,0])
    #s21right2.append(sp_right2[1,0])
    #s21.append(sp[1,0])

#plt.plot(f/1e9, 20*np.log10(np.abs(s21))) 
plt.plot(f/1e9, 20*np.log10(np.abs(s21left))) 
#plt.plot(f/1e9, 20*np.log10(np.abs(s21right)))
#plt.plot(f/1e9, 20*np.log10(np.abs(s21right2))) 
plt.xlabel("Frequency(GHz)")
plt.ylabel("S21(dB)")


# In[1106_InAs transmission line check]

path = 'C:/Users/20141/OneDrive - purdue.edu/PhD_research/resonators/reson1106_InAs/'
sp100nm = '100nm_s12.csv'
sp10nm = '10nm_s12.csv'
data = [sp100nm, sp10nm]
lab = ['Thick 100nm', 'Thick 10nm']

for i in range(len(data)):
    sp = data[i]
    dat = np.loadtxt(path+sp, delimiter=';', skiprows=1)
    plt.figure(figsize=(5,5))
    plt.plot(dat[:,0]/1e9, 20*np.log10(dat[:,1]),'.-', label=lab[i])
    plt.xlabel('Frequency(GHz)')
    plt.ylabel('S12(dB)')
    
# In[1109_ model circuit check]

lenlist = [5.17e-3, 10.34e-3, 100e-6]
Llist = [4.787e-7, 4.787e-7, 4.787e-4]
Clist = [1.951e-10, 1.951e-10, 1.951e-10]



# In[1109_ model circuit check]

lenlist = [5.17e-3, 10.34e-3]
Llist = [4.787e-7, 4.787e-7]
Clist = [1.951e-10, 1.951e-10]

# In[1109_ model circuit check]
for i in range(10):
    print(vph(l=lenlist, L=Llist, C=Clist, QWR=2*i+1)[2]/1e9)
    
# In[1109_ model circuit check]
path = 'C:/Users/20141/OneDrive - purdue.edu/PhD_research/resonators/reson1106_InAs/'
slist = ['S21_F4G_T5G_v2.csv', 'S21_F1G_T2G.csv', 'S21_F8G_T9G.csv', 'S21_F15G_T16G.csv']
flist = [np.linspace(5.002, 5.00204, 10000),np.linspace(1.66, 1.8, 10000),np.linspace(8.34, 8.35, 10000),np.linspace(15.005, 15.007, 10000)]
for i in range(len(slist)):
    plt.figure(figsize=(7,3))
    raw = np.loadtxt(path+slist[i], skiprows=1, delimiter=';')
    plt.plot(flist[i],raw[:,1],'.')
    plt.xlabel('Frequency(GHZ)')
    plt.ylabel('S21(dB)')    
    
    
    