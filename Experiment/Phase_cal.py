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
soccfg = xilinx.socconfig()[1]  # SoC configuration
soc = xilinx.socconfig()[0]   # SoC 

out_chs = [1]  # gen ch=1., ro_ch=1

# In[2]: Frequency range to be tested & helper functions

# Output frequency f0_v.
###################
# Try it yourself ! Change the output frequency.
###################

f0_start=4900
f0_step=1
expts=300
# expts=1

f0_v = np.arange(0,expts)*f0_step+f0_start

f0_v = soccfg.adcfreq(f0_v, gen_ch=1, ro_ch=1)

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



# In[5]: Averager Program
class SingleFreqProgram(AveragerProgram):
    def __init__(self,soccfg, cfg):
        super().__init__(soccfg, cfg)

    def initialize(self):
        cfg=self.cfg   
        
        # configure the readout lengths and downconversion frequencies
        for ch in range(2):
            self.declare_readout(ch=ch, length=self.cfg["readout_length"],
                                 freq=self.cfg["pulse_freq"])

        idata = 30000*np.ones(16*cfg["length"])

        for ch in self.cfg['out_ch']:
            self.declare_gen(ch=ch, nqz=1)
            self.add_pulse(ch=ch, name="measure", idata=idata)
        
        freq=soccfg.freq2reg(cfg["pulse_freq"])  # convert frequency to dac frequency
        self.trigger(pins=[0], t=0) # send a pulse on pmod0_0, for scope trigger
        for ch in self.cfg['out_ch']:
            self.set_pulse_registers(ch=ch, style="arb", freq=freq, phase=cfg["res_phase"], gain=cfg["pulse_gain"], 
                                     waveform="measure", mode="periodic")

        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        self.trigger(adcs=[0,1],adc_trig_offset=self.cfg["adc_trig_offset"])  # trigger the adc acquisition
        for ch in self.cfg['out_ch']:
            self.pulse(ch=ch, t=0) # play readout pulse
        self.wait_all() # control should wait until the readout is over
        self.sync_all(200)  # wait for measurement to complete

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

# In[6]: Decimated r/o

res=[]
config['reps'] = 1


# change this from 1 to 50
config['soft_avgs'] = 50


f0_start=5000
f0_step=0.001
expts=30

# =============================================================================
# f0_step=0.001
# expts=100
# f0_start=7437.54254 - f0_step*expts*0.5
# # expts=1 
# =============================================================================

f0_v = np.arange(0,expts)*f0_step+f0_start

f0_v = soccfg.adcfreq(f0_v, gen_ch=1, ro_ch=1)

# for f0 in [101]:
for f0 in f0_v:
    config['pulse_freq'] = f0
    prog =SingleFreqProgram(soccfg, config)
#     print(prog)
    data = prog.acquire_decimated(soc, load_pulses=True, progress=False, debug=False)
#     print(data)
    res.append([calculate_phase(d) for d in data])
    print_results(res[-1], f0)
fi_v=np.array([[a[0] for a in r] for r in res]).T

soc.reset_gens()

# In[7]: accumulated r/o
# change this from 10 to 1000
from scipy.optimize import least_squares

def phase_residuals(data,prediction):
    r = np.remainder(data-prediction+180,360)-180
    return r
    
def phase_model(x, f0):
    return np.remainder(x[0] - 360*x[1]*(f0), 360)

def phase_func(x, arg):
    resid = phase_residuals(arg, phase_model(x, f0_v))
    return resid

def phasecal(f0_start, f0_step, expts):
    
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
    
    config['reps'] = 1000
    
    config['soft_avgs'] = 1
    
    f0_v = np.arange(0,expts)*f0_step+f0_start
    f0_v = soccfg.adcfreq(f0_v, gen_ch=1, ro_ch=1)
    
    res = []
    adc1_ph = []
    
    for f0 in f0_v:
        config['pulse_freq'] = f0
        prog =SingleFreqProgram(soccfg, config)
        avg_data = prog.acquire(soc, load_pulses=True, progress=False, debug=False)
        data = [[prog.di_buf[i]/config['readout_length'], prog.dq_buf[i]/config['readout_length']] for i in range(2)]
        res.append([calculate_phase(d) for d in data])
        print_results(res[-1], f0)
        adc1_ph.append(res[-1][1][0])
        
    fi_v=np.array([[a[0] for a in r] for r in res]).T
    soc.reset_gens()
    

    fig, axs = plt.subplots(2,2)
    for ch in range(2):
    
        slopes = -1*(fi_v[ch,1:]-fi_v[ch,:-1])/(360*(f0_v[1:]-f0_v[:-1]))
    #     print(slopes)
        x0 = np.zeros(2)
        x0[1] = np.median(slopes)
        x0[0] = np.remainder(np.median(phase_residuals(fi_v[ch],phase_model(x0, f0_v))),360)
        print("initial estimate: %.2f deg shift, %.2f us delay"% tuple(x0))
    
        fit = least_squares(phase_func, x0, args=(fi_v[ch],))
    #     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale='jac')
    #     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale=(1,1e-5))
        fit.x[0] = np.remainder(fit.x[0],360)
    #     print(fit.status)
        print("after minimization: %.2f deg shift, %.2f us delay"% tuple(fit.x))
    
        plot = axs[0,ch]
        plot.set_title(r"$\phi$ vs $f$")
        plot.set_ylabel(r"$\phi$ (degrees)")
        plot.set_xlabel(r"$f$ (MHz)")
        plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
        plot.plot(f0_v, phase_model(fit.x, f0_v))
    
        plot = axs[1,ch]
        plot.set_ylabel(r"residual (degrees)")
        plot.set_xlabel(r"$f$ (MHz)")
        
        plot.plot(f0_v, phase_func(fit.x, fi_v[ch]), marker='.', linestyle="None",color="Red")
        
    return adc1_ph
# In[]

vna_ph = phasecal(4700, 1, 600)

# In[]

ph2 = phasecal(5030, 0.0025, 300)


# In[8] plot
fig, axs = plt.subplots(2,1)

for ch in range(2):
    plot = axs[ch]
    plot.plot(f0_v,fi_v[ch])
    plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
    plot.set_title(r"$\phi$ vs $f$ : ADC_{}".format(ch))
    plot.set_ylabel(r"$\phi$ (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
plt.tight_layout()   
# In[9] fit

def phase_residuals(data,prediction):
    r = np.remainder(data-prediction+180,360)-180
    return r
    
def phase_model(x, f0):
    return np.remainder(x[0] - 360*x[1]*(f0), 360)

def phase_func(x, arg):
    resid = phase_residuals(arg, phase_model(x, f0_v))
    return resid

from scipy.optimize import least_squares
fig, axs = plt.subplots(2,2)
for ch in range(2):

    slopes = -1*(fi_v[ch,1:]-fi_v[ch,:-1])/(360*(f0_v[1:]-f0_v[:-1]))
#     print(slopes)
    x0 = np.zeros(2)
    x0[1] = np.median(slopes)
    x0[0] = np.remainder(np.median(phase_residuals(fi_v[ch],phase_model(x0, f0_v))),360)
    print("initial estimate: %.2f deg shift, %.2f us delay"% tuple(x0))

    fit = least_squares(phase_func, x0, args=(fi_v[ch],))
#     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale='jac')
#     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale=(1,1e-5))
    fit.x[0] = np.remainder(fit.x[0],360)
#     print(fit.status)
    print("after minimization: %.2f deg shift, %.2f us delay"% tuple(fit.x))

    plot = axs[0,ch]
    plot.set_title(r"$\phi$ vs $f$")
    plot.set_ylabel(r"$\phi$ (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
    plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
    plot.plot(f0_v, phase_model(fit.x, f0_v))

    plot = axs[1,ch]
    plot.set_ylabel(r"residual (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
    
    plot.plot(f0_v, phase_func(fit.x, fi_v[ch]), marker='.', linestyle="None",color="Red")
#     plot.plot(f0_v, phase_func(x0, (fi_v[ch])), marker='.', linestyle="None",color="Red")


'''

'''



# =============================================================================
# 
# # In[3] Measurement methods
# 
# def measure_phase_decimated(soc, f0, fi0, out_array):
#     """
#     Measure the phase shift for a given frequency using decimated readout, and write it into the provided array.
#     Assumes that the tProc is already loaded with a program that reads the frequency and phase from memory locations 123 and 124.
#     """
#     readout_length = 1000
# 
#     f0_dac_int = soccfg.freq2reg(f=f0) # Frequency of signal generator
#     
#     soc.tproc.single_write(addr=123, data=f0_dac_int)
# 
#     fi0_int = soccfg.deg2reg(fi0, gen_ch=1)
#     soc.tproc.single_write(addr=124, data=fi0_int)
# 
#     # Configure readout.
#     for ch in range(2):
#         soc.configure_readout(ch=ch, output="product", frequency=f0)
# #         soc.config_avg(ch=ch, address=0, length=readout_length, enable=True)
#         soc.config_buf(ch=ch, address=0, length=readout_length, enable=True)
#     
#     # Start tProc.
#     soc.tproc.start()
# 
#     time.sleep(0.1)
#     
#     # Get data.
#     data = [soc.get_decimated(ch=ch, length=readout_length) for ch in range(2)]
# #     print(data)
# 
#     # Process data.
#     res = [calculate_phase(d) for d in data]
# 
#     # Save results.
#     out_array[:] = [a[0] for a in res]
#     print_results(res, f0)
#     
# # In[4]
# 
# nsamp = 10
# 
# prog = QickProgram(soccfg)
# prog.memri(1,1,123,"freq")
# prog.memri(1,2,124,"phase")
# prog.regwi(1,3,32000,"gain")
# prog.regwi(1,4,nsamp,"nsamp: generator will consume 16*nsamp DAC values.")
# prog.regwi(1,5,0x4,"b00100 -> phrst = 0, stdysel = 0, mode = 1, outsel = 00")
# prog.bitwi(1,5,5,"<<",16)
# prog.bitw(1,4,4,"|",5)
# prog.regwi(1,6,0,"start time")
# prog.synci(1000)
# prog.comment("Set trigger.")
# prog.regwi(0,1,0x1,"Send a pulse on pmod 0_0 (pin 1 on J48 on the ZCU111).")
# prog.seti(0,0,1,0,"Start the pulse.")
# prog.seti(0,0,0,100,"End the pulse after 100 ticks (260 ns).")
# 
# prog.comment("Program signal generators in periodic mode.") #TODO: make comments sppear in the printout?
# prog.comment("Channels 4, 5, 6, 7 -> DAC 229 CH0/1/2/3.")
# for ch in out_chs:
#     # the channel number here is 1-indexed, since it's the tProc channel number for the gen
#     prog.set(ch+1,1,1,2,0,3,4,6)
# prog.synci(1000)
# prog.comment("Set trigger.")
# prog.regwi(0,1,0xc000,"Trigger both buffers.")
# prog.seti(0,0,1,0,"Start the pulse.")
# prog.seti(0,0,0,100,"End the pulse after 100 ticks (260 ns).")
# prog.end()
# 
# print(prog)
# 
# prog.load_program(soc)
# 
# # Constant, real envelope.
# # The length must be at least 16 times the nsamp parameter passed to the signal generators.
# xg_i = 30000*np.ones(16*nsamp)
# xg_q = np.zeros_like(xg_i)
# 
# for ch in out_chs:
#     # Set the DAC channels to be in 1st Nyquist zone mode
#     soc.set_nyquist(ch=ch,nqz=1)
#     # Upload waveform.
#     soc.load_pulse_data(ch=ch,idata=xg_i, qdata=xg_q,addr=0)
# 
# # Output phase.
# fi0 = 0
# 
# fi_v = np.zeros((2,len(f0_v)))
# 
# # Load program.
# for ii, f0 in enumerate(f0_v):
#     measure_phase_decimated(soc, f0, fi0, fi_v[:,ii])
# 
# soc.reset_gens()
# 
# '''
# doesn't work!
# '''
# =============================================================================

# In[2023 0414]:
    
f0_start =4000
f0_step=0.0005
expts=50
# expts=1

f0_v = np.arange(0,expts)*f0_step+f0_start

f0_v = soccfg.adcfreq(f0_v, gen_ch=1, ro_ch=0)

# In[decimated r/o]:
    
res=[]
config['reps'] = 1

# change this from 1 to 50
config['soft_avgs'] = 50

# for f0 in [101]:
for f0 in f0_v:
    config['pulse_freq'] = f0
    prog =SingleFreqProgram(soccfg, config)
#     print(prog)
    data = prog.acquire_decimated(soc, load_pulses=True, progress=False, debug=False)
#     print(data)
    res.append([calculate_phase(d) for d in data])
    print_results(res[-1], f0)
fi_v=np.array([[a[0] for a in r] for r in res]).T

soc.reset_gens()

# In[accumulated r/o]:
    
# change this from 10 to 1000
config['reps'] = 1000

config['soft_avgs'] = 1
res=[]
for f0 in f0_v:
    config['pulse_freq'] = f0
    prog =SingleFreqProgram(soccfg, config)
    avg_data = prog.acquire(soc, load_pulses=True, progress=False, debug=False)
    data = [[prog.di_buf[i]/config['readout_length'], prog.dq_buf[i]/config['readout_length']] for i in range(2)]
    res.append([calculate_phase(d) for d in data])
    print_results(res[-1], f0)
fi_v=np.array([[a[0] for a in r] for r in res]).T
soc.reset_gens()



# In[]

# Plot results.
# plt.figure(2)
fig, axs = plt.subplots(2,1)

for ch in range(2):
    plot = axs[ch]
    plot.plot(f0_v,fi_v[ch])
    plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
    plot.set_title(r"$\phi$ vs $f$")
    plot.set_ylabel(r"$\phi$ (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
    
# In[fit]:
    
def phase_residuals(data,prediction):
    r = np.remainder(data-prediction+180,360)-180
    return r
    
def phase_model(x, f0):
    return np.remainder(x[0] - 360*x[1]*(f0), 360)

def phase_func(x, arg):
    resid = phase_residuals(arg, phase_model(x, f0_v))
    return resid

from scipy.optimize import least_squares
fig, axs = plt.subplots(2,2)
for ch in range(2):

    slopes = -1*(fi_v[ch,1:]-fi_v[ch,:-1])/(360*(f0_v[1:]-f0_v[:-1]))
#     print(slopes)
    x0 = np.zeros(2)
    x0[1] = np.median(slopes)
    x0[0] = np.remainder(np.median(phase_residuals(fi_v[ch],phase_model(x0, f0_v))),360)
    print("initial estimate: %.2f deg shift, %.2f us delay"% tuple(x0))

    fit = least_squares(phase_func, x0, args=(fi_v[ch],))
#     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale='jac')
#     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale=(1,1e-5))
    fit.x[0] = np.remainder(fit.x[0],360)
#     print(fit.status)
    print("after minimization: %.2f deg shift, %.2f us delay"% tuple(fit.x))

    plot = axs[0,ch]
    plot.set_title(r"$\phi$ vs $f$")
    plot.set_ylabel(r"$\phi$ (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
    plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
    plot.plot(f0_v, phase_model(fit.x, f0_v))

    plot = axs[1,ch]
    plot.set_ylabel(r"residual (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
    
    plot.plot(f0_v, phase_func(fit.x, fi_v[ch]), marker='.', linestyle="None",color="Red")
#     plot.plot(f0_v, phase_func(x0, (fi_v[ch])), marker='.', linestyle="None",color="Red")

# In[accumulated r/o with delay compensation]:
    
# change this from 10 to 1000
config['reps'] = 1000

config['soft_avgs'] = 1
res=[]
for f0 in f0_v:
    config['pulse_freq'] = f0
    ph_delay = f0*41.92*360
    config["res_phase"]=soccfg.deg2reg(ph_delay)
    prog =SingleFreqProgram(soccfg, config)
    avg_data = prog.acquire(soc, load_pulses=True, progress=False, debug=False)
    data = [[prog.di_buf[i]/config['readout_length'], prog.dq_buf[i]/config['readout_length']] for i in range(2)]
    res.append([calculate_phase(d) for d in data])
    print_results(res[-1], f0)
fi_v=np.array([[a[0] for a in r] for r in res]).T
soc.reset_gens()


# In[]

# Plot results.
# plt.figure(2)
fig, axs = plt.subplots(2,1)

for ch in range(2):
    plot = axs[ch]
    plot.plot(f0_v,fi_v[ch])
    plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
    plot.set_title(r"$\phi$ vs $f$")
    plot.set_ylabel(r"$\phi$ (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
    
for ch in range(2):

    slopes = -1*(fi_v[ch,1:]-fi_v[ch,:-1])/(360*(f0_v[1:]-f0_v[:-1]))
#     print(slopes)
    x0 = np.zeros(2)
    x0[1] = np.median(slopes)
    x0[0] = np.remainder(np.median(phase_residuals(fi_v[ch],phase_model(x0, f0_v))),360)
    print("initial estimate: %.2f deg shift, %.2f us delay"% tuple(x0))

    fit = least_squares(phase_func, x0, args=(fi_v[ch],))
#     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale='jac')
#     fit = least_squares(phase_func, x0, args=(fi_v[ch]), method='lm', x_scale=(1,1e-5))
    fit.x[0] = np.remainder(fit.x[0],360)
#     print(fit.status)
    print("after minimization: %.2f deg shift, %.2f us delay"% tuple(fit.x))

    plot = axs[0,ch]
    plot.set_title(r"$\phi$ vs $f$")
    plot.set_ylabel(r"$\phi$ (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
    plot.plot(f0_v,fi_v[ch], marker='.', linestyle="None",color="Red")
    plot.plot(f0_v, phase_model(fit.x, f0_v))

    plot = axs[1,ch]
    plot.set_ylabel(r"residual (degrees)")
    plot.set_xlabel(r"$f$ (MHz)")
    
    plot.plot(f0_v, phase_func(fit.x, fi_v[ch]), marker='.', linestyle="None",color="Red")
#     plot.plot(f0_v, phase_func(x0, (fi_v[ch])), marker='.', linestyle="None",color="Red")
