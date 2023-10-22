# In[0] Phase cal.
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('Z:/Mingi Kim/python_codes')
sys.path.append('Z:/Mingi Kim/FPGA/qick/qick_lib')
from qick import *
import LAB_v0 
import os
import time

# In[1] RFSoC
xilinx = LAB_v0.RFSoC()

# In[1.5]:
    
'''
Phase calibration for gen ch=1, ro_ch = 0
'''

soccfg = xilinx.socconfig()[1]  # SoC configuration
soc = xilinx.socconfig()[0]   # SoC 

out_chs = [1]  # gen ch=1., ro_ch=0

def calculate_phase(d):
    [xi,xq] = d
    x = xi +1j*xq

    # Average to improve calibration.
    xavg = np.mean(x)

    # Calculate calibration phase.
    fi = np.remainder(np.angle(xavg,deg=True)+360,360)
    return [fi, np.abs(xavg), np.std(x)]

def print_results(res, f0):
    print("freq_i = %f MHz, "%(f0) +
         "phi_i = (%.2f, %.2f) deg, " % tuple([res[i][0] for i in range(2)]) +
         "mag = (%.2f, %.2f), " % tuple([res[i][1] for i in range(2)]) +
         "RMS = (%.2f, %.2f) ADU" % tuple([res[i][2] for i in range(2)]))
#     print("freq_i = %f MHz, phi_i = (%.2f, %.2f) deg, mag = (%.2f, %.2f), RMS = (%.2f, %.2f) ADU" %(f0,*out_array,*A,*xrms))



# In[6]: Decimated r/o
config={"out_ch":out_chs,
        "reps":1, # --Fixed
        "res_phase":soccfg.deg2reg(0), # --Fixed
        
        "length":10, # [Clock ticks]
        # Try varying length from 10-100 clock ticks
        
        "readout_length":1000, # [Clock ticks]
        # Try varying readout_length from 50-1000 clock ticks

        "pulse_gain":32000, # [DAC units]
        # Try varying pulse_gain from 500 to 30000 DAC units

        "pulse_freq": 1400, # [MHz]
        # In this program the signal is up and downconverted digitally so you won't see any frequency
        # components in the I/Q traces below. But since the signal gain depends on frequency, 
        # if you lower pulse_freq you will see an increased gain.

        "adc_trig_offset": 1000, # [Clock ticks]
        # Try varying adc_trig_offset from 100 to 220 clock ticks

        "soft_avgs":1
        # Try varying soft_avgs from 1 to 200 averages

       }

res=[]
config['reps'] = 1


# change this from 1 to 50
config['soft_avgs'] = 50


f0_start=5000
f0_step=0.0007
expts=80

# =============================================================================
#f0_start=100
#f0_step=0.000250
#expts=40
# =============================================================================

f0_v = np.arange(0,expts)*f0_step+f0_start

f0_v = soccfg.adcfreq(f0_v, gen_ch=1, ro_ch=0)

# for f0 in [101]:
for f0 in f0_v:
    config['pulse_freq'] = f0
    prog =LAB_v0.SingleFreqProgram(soccfg, config)
#     print(prog)
    data = prog.acquire_decimated(soc, load_pulses=True, progress=False, debug=False)
#     print(data)
    res.append([calculate_phase(d) for d in data])
    print_results(res[-1], f0)
fi_v=np.array([[a[0] for a in r] for r in res]).T

soc.reset_gens()

# In[plot freq vs data]

plt.plot(f0_v, fi_v[0])
plt.plot(f0_v, fi_v[0], marker='.', linestyle="None", color='Red')
plt.title(r"$\phi$ vs $f$ : ADC_{}".format(0))
plt.ylabel(r"$\phi$ (degrees)")
plt.xlabel(r"$f$ (MHz)")

# In[9] fit the phase difference

def phase_residuals(data,prediction):
    r = np.remainder(data-prediction+180,360)-180
    return r
    
def phase_model(x, f0):
    return np.remainder(x[0] - 360*x[1]*(f0), 360)

def phase_func(x, arg):
    resid = phase_residuals(arg, phase_model(x, f0_v))
    return resid

# In[10]:
    
from scipy.optimize import least_squares
#fig, axs = plt.subplots(2,1)

ch = 0

slopes = -1*(fi_v[ch,1:]-fi_v[ch,:-1])/(360*(f0_v[1:]-f0_v[:-1]))
print(slopes)


x0 = np.zeros(2)
x0[1] = np.median(slopes)
x0[0] = np.remainder(np.median(phase_residuals(fi_v[ch],phase_model(x0, f0_v))),360)
print("initial estimate: %.2f deg shift, %.2f us delay"% tuple(x0))

# In[11]:
fit = least_squares(phase_func, x0, max_nfev=500, args=(fi_v[ch],))
#     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale='jac')
#     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale=(1,1e-5))
fit.x[0] = np.remainder(fit.x[0], 360)
print(fit.status)
print("\nafter minimization: %.2f deg shift, %.2f us delay"% tuple(fit.x))

# In[]
fig, axs = plt.subplots(2,1)

plot = axs[0]
plot.set_title(r"$\phi$ vs $f$")
plot.set_ylabel(r"$\phi$ (degrees)")
plot.set_xlabel(r"$f$ (MHz)")
plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
plot.plot(f0_v, phase_model(fit.x, f0_v))

plot = axs[1]
plot.set_ylabel(r"residual (degrees)")
plot.set_xlabel(r"$f$ (MHz)")

plot.plot(f0_v, phase_func(fit.x, fi_v[ch]), marker='.', linestyle="None",color="Red")
#     plot.plot(f0_v, phase_func(x0, (fi_v[ch])), marker='.', linestyle="None",color="Red")



# In[5]:
'''
After phase calibration : check for the phase calibrated data
'''
    
    
# In[6]: frequency setup

res=[]
config['reps'] = 1


# change this from 1 to 50
config['soft_avgs'] = 50


f0_start=5000
f0_step=3
expts=100

# =============================================================================
#f0_start=100
#f0_step=0.000250
#expts=40
# =============================================================================

f0_v = np.arange(0,expts)*f0_step+f0_start

f0_v = soccfg.adcfreq(f0_v, gen_ch=1, ro_ch=1)

# In[6]: Decimated r/o

# for f0 in [101]:
for f0 in f0_v:
    config['pulse_freq'] = f0
    prog =LAB_v0.SingleFreqProgram(soccfg, config)
#     print(prog)
    data = prog.acquire_decimated(soc, load_pulses=True, progress=False, debug=False)
#     print(data)
    res.append([calculate_phase(d) for d in data])
    print_results(res[-1], f0)
fi_v=np.array([[a[0] for a in r] for r in res]).T

soc.reset_gens()

# In[accumulated r/o with delay compensation]:
    
# change this from 10 to 1000
config['reps'] = 1000

config['soft_avgs'] = 1
res=[]
for f0 in f0_v:
    config['pulse_freq'] = f0

    prog =LAB_v0.SingleFreqProgram(soccfg, config)
    avg_data = prog.acquire(soc, load_pulses=True, progress=False, debug=False)
    data = [[prog.di_buf[i]/config['readout_length'], prog.dq_buf[i]/config['readout_length']] for i in range(2)]
    res.append([calculate_phase(d) for d in data])
    print_results(res[-1], f0)
fi_v=np.array([[a[0] for a in r] for r in res]).T
soc.reset_gens()


# In[]
fig, axs = plt.subplots(2,1)

plot = axs[0]
plot.set_title(r"$\phi$ vs $f$")
plot.set_ylabel(r"$\phi$ (degrees)")
plot.set_xlabel(r"$f$ (MHz)")
plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
plot.plot(f0_v, phase_model(fit.x, f0_v))

plot = axs[1]
plot.set_ylabel(r"residual (degrees)")
plot.set_xlabel(r"$f$ (MHz)")

plot.plot(f0_v, phase_func(fit.x, fi_v[ch]), marker='.', linestyle="None",color="Red")
#     plot.plot(f0_v, phase_func(x0, (fi_v[ch])), marker='.', linestyle="None",color="Red")
