# In[0]:
    
import numpy as np
import pandas as pd
import sys
import os
import time
import matplotlib.pyplot as plt

sys.path.append('Z:/general/LRlabcode/LRlab/Experiment')
#sys.path.append('Z:/general/LRlabcode/LRlab/Experiment/Instruments')
sys.path.append('Z:/Mingi Kim/python_codes/Experiment')


# In[1]:
    
import LAB_v0
from data_analysis_codes import *
from Measurement import visualS21, frequency_sweep, VNA, bfieldsweep
from fit_codes import *

from qick import *

# In[2]:
    
xilinx = LAB_v0.RFSoC()
soccfg = xilinx.socconfig()[1]  # SoC configuration
soc = xilinx.socconfig()[0]   # SoC 

# In[3]

'''
channel config: 
    
    DAC_A [1]----coupler(in)--------ADC_D [0]
    DAC_B [0]----coupler(-20dB)-----ADC_D [0]

'''

hw_cfg={
        "drive_ch":1,
        "probe_ch":0
        }

readout_cfg={
    "readout_length":soccfg.us2cycles(3.0, gen_ch=0),
    "f_probe":99.775+0.18, #MHz
    "res_phase":0,
    "adc_trig_offset": 130,
    "res_gain":10000
    }

drive_cfg={
    "sigma":soccfg.us2cycles(0.025, gen_ch=1),
    "pi_gain":11500,
    "pi2_gain":11500//2,
    "f_ge":4000,
    "relax_delay":500
    }




# In[helper functions]:
    
#helper functions
def hist(data=None, plot=True, ran=1.0):
    
    ig = data[0]
    qg = data[1]
    ie = data[2]
    qe = data[3]

    numbins = 200
    
    xg, yg = np.median(ig), np.median(qg)
    xe, ye = np.median(ie), np.median(qe)

    if plot==True:
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
        fig.tight_layout()

        axs[0].scatter(ig, qg, label='g', color='b', marker='*')
        axs[0].scatter(ie, qe, label='e', color='r', marker='*')
        axs[0].scatter(xg, yg, color='k', marker='o')
        axs[0].scatter(xe, ye, color='k', marker='o')
        axs[0].set_xlabel('I (a.u.)')
        axs[0].set_ylabel('Q (a.u.)')
        axs[0].legend(loc='upper right')
        axs[0].set_title('Unrotated')
        axs[0].axis('equal')
    """Compute the rotation angle"""
    theta = -np.arctan2((ye-yg),(xe-xg))
    """Rotate the IQ data"""
    ig_new = ig*np.cos(theta) - qg*np.sin(theta)
    qg_new = ig*np.sin(theta) + qg*np.cos(theta) 
    ie_new = ie*np.cos(theta) - qe*np.sin(theta)
    qe_new = ie*np.sin(theta) + qe*np.cos(theta)
    
    """New means of each blob"""
    xg, yg = np.median(ig_new), np.median(qg_new)
    xe, ye = np.median(ie_new), np.median(qe_new)
    
    #print(xg, xe)
    
    xlims = [xg-ran, xg+ran]
    ylims = [yg-ran, yg+ran]

    if plot==True:
        axs[1].scatter(ig_new, qg_new, label='g', color='b', marker='*')
        axs[1].scatter(ie_new, qe_new, label='e', color='r', marker='*')
        axs[1].scatter(xg, yg, color='k', marker='o')
        axs[1].scatter(xe, ye, color='k', marker='o')    
        axs[1].set_xlabel('I (a.u.)')
        axs[1].legend(loc='lower right')
        axs[1].set_title('Rotated')
        axs[1].axis('equal')

        """X and Y ranges for histogram"""
        
        ng, binsg, pg = axs[2].hist(ig_new, bins=numbins, range = xlims, color='b', label='g', alpha=0.5)
        ne, binse, pe = axs[2].hist(ie_new, bins=numbins, range = xlims, color='r', label='e', alpha=0.5)
        axs[2].set_xlabel('I(a.u.)')       
        
    else:        
        ng, binsg = np.histogram(ig_new, bins=numbins, range = xlims)
        ne, binse = np.histogram(ie_new, bins=numbins, range = xlims)

    """Compute the fidelity using overlap of the histograms"""
    contrast = np.abs(((np.cumsum(ng) - np.cumsum(ne)) / (0.5*ng.sum() + 0.5*ne.sum())))
    tind=contrast.argmax()
    threshold=binsg[tind]
    fid = contrast[tind]
    axs[2].set_title(f"Fidelity = {fid*100:.2f}%")

    return fid, threshold, theta

# In[Time of flight measurement]

'''
TOF : measures the time at which the measurement pulse appears in the adc buffer.
we only want to start capturing data from this point in time onwards --> will be stored in readout_cfg["adc_trig_offset"]

'''


class LoopbackProgram(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   
#         self.declare_gen(ch=cfg["jpa_ch"], nqz=1) #JPA
        self.declare_gen(ch=cfg["probe_ch"], nqz=1) #Readout
#         self.declare_gen(ch=cfg["qubit_ch"], nqz=2) #Qubit
#         self.declare_gen(ch=cfg["storage_ch"], nqz=2) #Storage

        for ch in [0,1]: #configure the readout lengths and downconversion frequencies
            self.declare_readout(ch=ch, length=cfg["readout_length"],
                                 freq=cfg["frequency"], gen_ch=cfg["probe_ch"])

        freq=self.freq2reg(cfg["frequency"], gen_ch=cfg["probe_ch"], ro_ch=0)  # convert frequency to dac frequency (ensuring it is an available adc frequency)
        self.set_pulse_registers(ch=cfg["probe_ch"], style="const", freq=freq, phase=0, gain=cfg["pulse_gain"],
                                length=cfg["pulse_length"])
        
        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        cfg=self.cfg   
        self.measure(pulse_ch=cfg["probe_ch"], 
             adcs=[0],
             adc_trig_offset=cfg["adc_trig_offset"],
             t=0,
             wait=True,
             syncdelay=self.us2cycles(cfg["relax_delay"]))        

expt_config={
        "reps":1, # --Fixed
        "pulse_length":600, # [Clock ticks]
        "readout_length":1000, # [Clock ticks]
        "pulse_gain":5000, # [DAC units]
        "frequency": 100, # [MHz]
        "adc_trig_offset": 150, # [Clock ticks]
        "soft_avgs":1000
       }


# In[TOF measure]

config={**hw_cfg,**readout_cfg,**drive_cfg,**expt_config}
prog =LoopbackProgram(soccfg, config)
adc1, adc2 = prog.acquire_decimated(soc, load_pulses=True, progress=True, debug=False)

# In[ Plot results.]
plt.figure(figsize=(15,10))
plt.subplot(111, title=f"Averages = {config['soft_avgs']}", xlabel="Clock ticks", ylabel="Transmission (adc levels)")
plt.plot(adc1[0], label="I value; ADC 0")
plt.plot(adc1[1], label="Q value; ADC 0")
#plt.plot(adc2[0], label="I value; ADC 1")
#plt.plot(adc2[1], label="Q value; ADC 1")
plt.legend()
plt.axvline(readout_cfg["adc_trig_offset"])


# In[single tone spectroscopy]
'''
Measures the single shot readout fidelity of the system. We acquire single shot (I, Q) readout values 
by first preparing the qubit in its ground (blue dots) a certain number of times (e.g. 5000 shots)
and then preparing the qubit in its excited state (red) the same number of times.
We then extract two parameters which are used to optimize the associated readout fidelity:
    the roation angle of the IQ blobs and the threshold that classified the two qubit states (g & e)
We store these two parameters here readout_cfg["res_phase"] - angle and readout_cfg["threshold"]

'''


class SingleShotProgram(RAveragerProgram):
    def initialize(self):
        cfg=self.cfg
        
        self.declare_gen(ch=cfg["probe_ch"], nqz=1) #Readout
        self.declare_gen(ch=cfg["drive_ch"], nqz=2) #Qubit (Andreev states)
        for ch in [0,1]: #configure the readout lengths and downconversion frequencies
            self.declare_readout(ch=ch, length=cfg["readout_length"],
                                 freq=cfg["f_probe"], gen_ch=cfg["probe_ch"])

        cfg["start"]=0
        cfg["step"]=cfg["pi_gain"]
        cfg["reps"]=cfg["shots"]
        cfg["expts"]=2
        
        self.q_rp=self.ch_page(self.cfg["drive_ch"])     # get register page for qubit_ch
        self.r_gain=self.sreg(cfg["drive_ch"], "gain")   # get frequency register for qubit_ch    
        
        f_res=self.freq2reg(cfg["f_probe"], gen_ch=cfg["probe_ch"], ro_ch=0) # conver f_res to dac register value
        f_ge=self.freq2reg(cfg["f_ge"], gen_ch=cfg["drive_ch"])
        
        # add qubit and readout pulses to respective channels
        self.add_gauss(ch=cfg["drive_ch"], name="drive", sigma=cfg["sigma"], length=cfg["sigma"]*4)
        self.set_pulse_registers(ch=cfg["drive_ch"], style="arb", freq=f_ge, phase=0, gain=cfg["start"], 
                                 waveform="drive")
        self.set_pulse_registers(ch=cfg["probe_ch"], style="const", freq=f_res, phase=cfg["res_phase"], gain=cfg["res_gain"], 
                                 length=cfg["readout_length"])

        self.sync_all(self.us2cycles(500))
    
    def body(self):
        self.pulse(ch=self.cfg["drive_ch"])  #play probe pulse
        self.sync_all(self.us2cycles(0.05)) # align channels and wait 50ns

        #trigger measurement, play measurement pulse, wait for qubit to relax
        self.measure(pulse_ch=self.cfg["probe_ch"], 
             adcs=[0],
             adc_trig_offset=self.cfg["adc_trig_offset"],
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"]))
    
    def update(self):
        self.mathi(self.q_rp, self.r_gain, self.r_gain, '+', self.cfg["step"]) # update frequency list index
        
    def acquire(self,soc, load_pulses=True, progress=False, debug=False):
        super().acquire(soc, load_pulses=load_pulses, progress=progress, debug=debug)
        return self.collect_shots()
        
    def collect_shots(self):
        shots_i0=self.di_buf[0].reshape((self.cfg["expts"],self.cfg["reps"]))/self.cfg['readout_length']
        shots_q0=self.dq_buf[0].reshape((self.cfg["expts"],self.cfg["reps"]))/self.cfg['readout_length']
        shots_i1=self.di_buf[1].reshape((self.cfg["expts"],self.cfg["reps"]))/self.cfg['readout_length']
        shots_q1=self.dq_buf[1].reshape((self.cfg["expts"],self.cfg["reps"]))/self.cfg['readout_length']
        return shots_i0,shots_q0,shots_i1,shots_q1
        
    def analyze(self, shots_i, shots_q):
        plt.subplot(111, xlabel='I', ylabel='Q', title='Single Shot Histogram')
        plt.plot(shots_i[0],shots_q[0],'.',label='g')
        plt.plot(shots_i[1],shots_q[1],'.',label='e')
        plt.legend()
        plt.gca().set_aspect('equal', 'datalim')

# In[experiment]

expt_cfg={
        "shots":5000, "res_phase":0
       }
config={**hw_cfg,**readout_cfg,**drive_cfg,**expt_cfg} #combine configs

ssp=SingleShotProgram(soccfg, config)
di0, dq0, di1, dq1 = ssp.acquire(soc, load_pulses=True,progress=True, debug=False)


# In[plot]
fid, threshold, angle = hist(data=[di0[0], dq0[0], di0[1], dq0[1]],  plot=True, ran=600)
#print('Optimal fidelity after rotation = %.3f' % fid)
readout_cfg["res_phase"]=soccfg.deg2reg(-angle*180/pi, gen_ch=0)
readout_cfg["threshold"]=10 #round(threshold)

# In[Pulse probe spectroscopy]

'''
IN the demo, it measures the qubit frequency f_ge

RAverageProgram --> sweeps a parameter directly on the processor rather than in a Python loop.
Because the whole sweep is done on the processor, there is less downtime (especially for fast experiments)

'''

class PulseProbeSpectroscopyProgram(RAveragerProgram):
    def initialize(self):
        cfg=self.cfg
        
        self.declare_gen(ch=cfg["probe_ch"], nqz=1) #Readout (probe)
        self.declare_gen(ch=cfg["drive_ch"], nqz=2) #ABS (drive)
        for ch in [0,1]: #configure the readout lengths and downconversion frequencies
            self.declare_readout(ch=ch, length=cfg["readout_length"],
                                 freq=cfg["f_probe"], gen_ch=cfg["probe_ch"])

        self.q_rp=self.ch_page(self.cfg["drive_ch"])     # get register page for drive_ch
        self.r_freq=self.sreg(cfg["drive_ch"], "freq")   # get frequency register for drive_ch    
        
        f_res=self.freq2reg(cfg["f_probe"], gen_ch=cfg["probe_ch"], ro_ch=0) # conver f_res to dac register value

        self.f_start =self.freq2reg(cfg["start"], gen_ch=cfg["drive_ch"])  # get start/step frequencies
        self.f_step =self.freq2reg(cfg["step"], gen_ch=cfg["drive_ch"])

        # add qubit and probe pulses to respective channels
        self.set_pulse_registers(ch=cfg["drive_ch"], style="const", freq=self.f_start, phase=0, gain=cfg["drive_gain"], 
                                 length=cfg["probe_length"])
        self.set_pulse_registers(ch=cfg["probe_ch"], style="const", freq=f_res, phase=cfg["res_phase"], gain=cfg["res_gain"], 
                                 length=cfg["readout_length"])
        
        self.sync_all(self.us2cycles(1))
    
    def body(self):
        self.pulse(ch=self.cfg["drive_ch"])  #play probe pulse
        self.sync_all(self.us2cycles(0.05)) # align channels and wait 50ns

        #trigger measurement, play measurement pulse, wait for qubit to relax
        self.measure(pulse_ch=self.cfg["probe_ch"], 
             adcs=[0],
             adc_trig_offset=self.cfg["adc_trig_offset"],
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"]))

    def update(self):
        self.mathi(self.q_rp, self.r_freq, self.r_freq, '+', self.f_step) # update fre

# In[]

expt_cfg={"start":4741, "step":0.01, "expts":400, "reps": 20,"rounds":50,
          "probe_length":soccfg.us2cycles(2.0, gen_ch=hw_cfg["drive_ch"]), "drive_gain":10
         }
config={**hw_cfg,**readout_cfg,**drive_cfg,**expt_cfg} #combine configs

qspec=PulseProbeSpectroscopyProgram(soccfg, config)
expt_pts, avgi, avgq = qspec.acquire(soc, threshold=readout_cfg["threshold"],load_pulses=True,progress=True, debug=False)

# In[]

plt.subplot(111,title="Qubit Spectroscopy", xlabel="Qubit Frequency (GHz)", ylabel="Qubit Population")
plt.plot(expt_pts, avgi[0][0],'o-')
#plt.axvline(drive_cfg["f_ge"]);

# In[Length rabi]

class LengthRabiProgram(AveragerProgram):
    def initialize(self):
        cfg=self.cfg

        self.declare_gen(ch=cfg["probe_ch"], nqz=1) #Readout
        self.declare_gen(ch=cfg["drive_ch"], nqz=2) #Qubit
        for ch in [0,1]: #configure the readout lengths and downconversion frequencies
            self.declare_readout(ch=ch, length=cfg["readout_length"],
                                 freq=cfg["f_probe"], gen_ch=cfg["probe_ch"])
        
        f_probe=self.freq2reg(cfg["f_probe"], gen_ch=cfg["probe_ch"], ro_ch=0) # conver f_res to dac register value
        f_ge=self.freq2reg(cfg["f_ge"], gen_ch=cfg["drive_ch"])

        # add qubit and readout pulses to respective channels
        self.set_pulse_registers(ch=cfg["drive_ch"], style="const", freq=f_ge, phase=0, gain=cfg["drive_gain"], 
                                 length=cfg["pulse_length"])
        self.set_pulse_registers(ch=cfg["probe_ch"], style="const", freq=f_res, phase=cfg["res_phase"], gain=cfg["res_gain"], 
                                 length=cfg["readout_length"])

        self.synci(200)
        
    def body(self):
        cfg=self.cfg
        self.pulse(ch=cfg["probe_ch"])  #play probe pulse
        self.sync_all(self.us2cycles(0.05)) # align channels and wait 50ns

        #trigger measurement, play measurement pulse, wait for qubit to relax
        self.measure(pulse_ch=cfg["probe_ch"], 
             adcs=[0],
             wait=True,
             syncdelay=self.us2cycles(cfg["relax_delay"]))

# In[expt]
expt_cfg={
       "qubit_gain":2000,
        "start":4, "step":1, "expts":200, "reps": 400,
       }
config={**hw_cfg,**readout_cfg,**qubit_cfg,**expt_cfg} #combine configs
expt_pts=[expt_cfg["start"] + ii*expt_cfg["step"] for ii in range(expt_cfg["expts"])]


# In[]

probe_list = np.round(np.linspace(2000, 2500, 10))
readout_list = np.round(np.linspace(1500, 3500, 10))
drive_list = np.round(np.linspace(2100,2300, 10 ))

for drive in drive_list:
    config["f_ge"] = drive
    
    for probe in probe_list:
        config["f_probe"]=probe
        
        rabi=LengthRabiProgram(soccfg, config)
        avgi,avgq = rabi.acquire(soc, threshold=readout_cfg["threshold"], load_pulses=True, progress=False,debug=False)
        results.append(avgi[0][0])


    
    