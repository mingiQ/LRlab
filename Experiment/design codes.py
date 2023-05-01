# In[0]:
import numpy as np
from numpy import cos, sin, tan, arctan, arcsin, arccos, pi
import matplotlib.pyplot as plt
from pyautocad import Autocad, APoint, aDouble
import sys
sys.path.append('C:/Users/20141/OneDrive - purdue.edu/PhD_research/resonators/python_codes/')
from resonator_params import *
import os

# In[1]: cpw design functions

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