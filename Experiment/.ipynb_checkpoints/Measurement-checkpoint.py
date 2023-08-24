# In[0]: necessary pkgs
import sys
import numpy as np
import time
import matplotlib.pyplot as plt

sys.path.append('Z:/Mingi Kim/python_codes')
import LAB_v0
import data_analysis_codes

# =============================================================================
# # In[1]: import xilinx SoC
#     
# xilinx = LAB_v0.RFSoC()
# 
# # In[1.5]: import RF instruments (anritsu, SoC and SoC Configuration)
# 
# vna = LAB_v0.MS2038('192.168.0.105')
# 
# soccfg = xilinx.socconfig()[1]  # SoC configuration
# soc = xilinx.socconfig()[0]   # SoC 
# =============================================================================

# In[measurement function]

def S21trans(start, stop, pnts, IF, avg, meas_t, path, name):
    '''
    Anritsu VNA code
    ----------
    start : start freq (GHz)
    stop : stop freq (GHz)
    IF : Intermediate Freq
    avg : # of avging
    meas_t : sleep time
    path : directory
    name : name

    Returns
    -------
    s2p files

    '''
    vna.avgcount(avg)                       # average count 
    vna.IFBW(IF)                         # IFBW 2kHz
    vna.start(start)
    vna.stop(stop)
    vna.points(pnts)
    
    filename = path+name
    
    vna.avgclear()
    time.sleep(meas_t)
    vna.save_s2p(filename)
    
def visualS21(start, stop, pnts, IF, avg, meas_t, path, name):
    S21trans(start, stop, pnts, IF, avg, meas_t, path, name)
    dat = data_extract(path, name)
    plt.figure(figsize=(8,4))
    plt.plot(dat[0], dat[1])
    plt.xlabel("Frequency(GHz)")
    plt.ylabel("S21(dB)")
    
def visualS21_plotall(start, stop, pnts, IF, avg, meas_t, path, name):
    
    S21trans(start, stop, pnts, IF, avg, meas_t, path, name)
    dat = data_extract(path, name)

    freq = dat[0]
    s21 = dat[1]
    phases = dat[2]
    
    
    plt.figure(figsize=(20,15))
    plt.subplot(221,title="DUT", xlabel="Frequency (GHz)", ylabel="S21(dB)")
    plt.plot(freq, s21)
    #max_freq=fpts[np.argmax(amps)]
    
    plt.subplot(222, title="phase", xlabel='Frequency', ylabel='Phase(degree)')
    plt.plot(freq, phases, '-')
    
    plt.subplot(212, projection='polar')
    #plt.axes(projection='polar')
    plt.plot(phases*np.pi/180, 10**(s21/20), '.')
    plt.grid(True)
    plt.show() 
    
def frequency_sweep(f_init, f_fin, scan_range, pnts, IF, avg, t, name, wdir):

    freq_list = np.arange(f_init, f_fin, scan_range/2)
    
    for freq in freq_list:
        f_start = (freq-scan_range/2)*1e9
        f_stop = (freq+scan_range/2)*1e9
    
    visualS21(start=f_start, stop=f_stop, pnts=pnts, IF=IF, avg=avg, meas_t=t, path=wdir, name=name)
    print("{}GHz scan completed\n".format(f_start))
    
    

def VNA(fstart, fstop, fstep, pgain, configs, wdir):
    
    
    '''
    RFSoC transmission measurement (ADC) ---> measure ADC level (arb)
    
    '''
    hw_cfg={"ro_ch":1,   # ADC_D
        "res_ch":1    # DAC_B
       }
    readout_cfg={
        "readout_length":soccfg.us2cycles(3.0, gen_ch=0), # [Clock ticks]
        "res_phase": 0,
        "adc_trig_offset": 275, # [Clock ticks]
        "res_gain_start":pgain
        }
    expt_cfg={"reps":500, "relax_delay":10,
              "start":fstart, "stop":fstop, "step":fstep,
              "gain_step" : 0.1
             }
    
    # reps: avgs, start: (MHz), step:(MHz), experiments (# of shots)
    
    config={**hw_cfg,**readout_cfg,**expt_cfg} #combine configs
    
    fpts=np.arange(expt_cfg["start"], expt_cfg["stop"], expt_cfg["step"])
    amps = []
    phases = []
    inphase = []
    quad = []
    
    print("measurement start!")
    time_init = time.time()
    for f in fpts:
        config["frequency"]=f
    
        rspec=LAB_v0.SingleToneSpectroscopyProgram(soccfg, config)
        avgi,avgq = rspec.acquire(soc, load_pulses=True)
        amp = np.abs(avgi[0][0]+1j*avgq[0][0])
        ph = np.remainder(np.angle(avgi[0][0]+1j*avgq[0][0],deg=True)+360,360)
        inphase.append(avgi[0][0])
        quad.append(avgq[0][0])
        amps.append(amp)
        phases.append(ph)
        
    meas_time = time.time()-time_init
    print("measurement done! ({}sec passed) ".format(meas_time))
    
    amps=np.array(amps)
    Is = np.array(inphase)
    Qs = np.array(quad)
    
    filename = '{}MHz to {}MHz pow {}adc_level'.format(fstart, fstop, pgain)+configs
    
    data = np.hstack((np.transpose([fpts]), np.transpose([amps]), np.transpose([phases]), np.transpose([Is]), np.transpose([Qs])))
    np.savetxt(wdir+filename, data, delimiter=',', header='freq,Amp,Phase,I,Q\nMHz,V,Deg,V,V', comments='')
    
        
    plt.figure(figsize=(20,15))
    plt.subplot(221,title="DUT", xlabel="Frequency (MHz)", ylabel="Amp. (adc level)")
    plt.plot(fpts, 20*np.log10(amps),'go-')
    max_freq=fpts[np.argmax(amps)]
    
    plt.subplot(222, title="phase", xlabel='Frequency', ylabel='Phase(degree)')
    plt.plot(fpts, phases, '-')
    
    plt.subplot(212, projection='polar')
    #plt.axes(projection='polar')
    plt.plot(phases, amps, '.')
    plt.grid(True)
    plt.show()

# =============================================================================
# # In[5] measurement code
# 
# hw_cfg={"ro_ch":0,   # ADC_D
#     "res_ch":1    # DAC_A
#    }
# readout_cfg={
#     "readout_length":soccfg.us2cycles(15.0, gen_ch=1), # [Clock ticks]
#     "res_phase": 0,
#     "adc_trig_offset": 275, # [Clock ticks]
#     "res_gain_start":3000
#     }
# expt_cfg={"reps":100, "relax_delay":15,
#           "start":4300, "step":0.1, "expts":500,
#           "center":4330, "span":100, "points":5000,
#           "gain_step" : 1, "soft_avg":1
#          }
# config={**hw_cfg,**readout_cfg,**expt_cfg} #combine configs
# 
# #fpts=expt_cfg["start"] + expt_cfg["step"]*np.arange(expt_cfg["expts"])
# fpts = np.linspace(expt_cfg["center"]-expt_cfg["span"]/2, expt_cfg["center"]+expt_cfg["span"]/2, expt_cfg["points"])
# gpts = config["res_gain_start"] + expt_cfg["gain_step"]*np.arange(expt_cfg["expts"])
# 
#     
#     
# def measure(center, span, plot_title, path):
#     expt_cfg["center"] = center
#     expt_cfg["span"] = span
#     fpts = np.linspace(expt_cfg["center"]-expt_cfg["span"]/2, expt_cfg["center"]+expt_cfg["span"]/2, expt_cfg["points"])
#     amps=[]
#     phases = []
#     
#     print("measurement start!")
#     time_init = time.time()
#     for f in fpts:
#         config["frequency"]=f
#     
#         rspec=LAB_v0.SingleToneSpectroscopyProgram(soccfg, config)
#         avgi,avgq = rspec.acquire(soc, load_pulses=True)
#         amp = np.abs(avgi[0][0]+1j*avgq[0][0])
#         ph = np.angle(avgi[0][0]+1j*avgq[0][0],deg=True)
#         amps.append(amp)
#         phases.append(ph)
#         
#     meas_time = time.time()-time_init
#     print("measurement done! ({}sec passed) ".format(meas_time))
#     
#     amps=np.array(amps)
#     phases = np.array(phases)
#     
#     spectrum_data = np.vstack((fpts, amps, phases))
#     np.savetxt(fname=path + plot_title+'csv', X=spectrum_data, delimiter=',', newline='\n')
#     
#     plt.figure(figsize=(20,15))
#     plt.subplot(221,title=plot_title, xlabel="Frequency (MHz)", ylabel="Amp. (adc level)")
#     plt.plot(fpts, amps,'go-')
#     
#     plt.subplot(222, title="phase", xlabel='Frequency', ylabel='Phase(degree)')
#     plt.plot(fpts, phases, '.-')
#     
#     plt.subplot(212, projection='polar')
#     plt.plot(phases*np.pi/180, amps, '.')
#     plt.grid(True)
#     plt.show()
#     print("resonance frequency: {}GHz".format(fpts[np.argmin(amps)]/1e3))
#     return fpts[np.argmin(amps)]
#     
# def power_sweep(center, span, plot_title, gain_start, gain_stop, gain_points):
#     gpts = np.linspace(gain_start, gain_stop, gain_points)
#     for g in gpts:
#         config["res_gain_start"] = g
#         measure(center, span, plot_title)
# 
# # In[1]: instruments
# 
# #vna = LAB_v0.MS2038('192.168.0.105')
# 
# =============================================================================

