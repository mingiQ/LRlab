# In[0]: necessary pkgs
import sys
import numpy as np
import time
import matplotlib.pyplot as plt

sys.path.append('Z:/general/LRlabcode/LRlab')
import Experiment.LAB_v0
from Experiment.data_analysis_codes import *

# In[measurement function]

########################################################################################
'''
Anritsu measurement

'''

def S21trans(vna, start, stop, pnts, IF, avg, meas_t, path, name):
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
    vna.sweep('single')
    time.sleep(meas_t)
    vna.save_s2p(filename)
    
def visualS21(vna, start, stop, pnts, IF, avg, meas_t, path, name):
    
    S21trans(vna, start, stop, pnts, IF, avg, meas_t, path, name)
    
    if vna != 'cmt':
        dat = data_extract(path, name, 'anri')
       
    elif vna == 'cmt':
        dat = data_extract(path, name, 'cmt')
   

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
    
    return freq, s21, phases
    
    
def frequency_sweep(vna, f_init, f_fin, scan_range, pnts, IF, avg, t, name, wdir):

    freq_list = np.arange(f_init, f_fin, scan_range/2)
    
    for freq in freq_list:
        f_start = (freq-scan_range/2)*1e9
        f_stop = (freq+scan_range/2)*1e9
    
        visualS21(vna=vna, start=f_start, stop=f_stop, pnts=pnts, IF=IF, avg=avg, meas_t=t, path=wdir, name=name+'near_{}GHz.s2p'.format(freq))
        print("{}Hz scan completed\n".format(f_start))
    
###############################################################################################

'''
RF SoC measurement

'''

def ADCcal(soc, config, offset, num):

    avglistv2 = []
    stdlistv2 = []
    attlistv2 = []
    
    for i in range(num):
        
        att.set_ch(1,i+offset)
        attlistv2.append(i+offset)
        #print(i+ 'att \n')
        
        prog=LAB_v0.LoopbackProgram_sendpulse_readpulse(soccfg, config)
        iq_list = prog.acquire_decimated(soc, load_pulses=True, progress=True, debug=False)
        
        mag_sq = iq_list[0][0][200-100:200+100]**2 + iq_list[0][1][200-100:200+100]**2
        mag = np.sqrt(mag_sq)
        
        avgmag = np.average(mag)
        stdmag = np.std(mag)
        
        print(i, avgmag, stdmag)
        
        avglistv2.append(avgmag)
        stdlistv2.append(stdmag)
    
    plt.figure(figsize=(10,5))
    plt.plot(attlistv2, avglistv2)
    plt.figure(figsize=(10,5))
    plt.plot(attlistv2, stdlistv2)
    plt.figure(figsize=(10,5))
    plt.plot(attlistv2, np.array(avglistv2)/np.array(stdlistv2))
    plt.title('SNR = avg/std')
    
    return (avglistv2, stdlistv2, attlistv2)
    

def VNA(soccfg, soc, gen, ro, rolen, fstart, fstop, fstep, pgain, configs, cal, wdir):
    
    '''
    RFSoC transmission measurement (ADC) ---> measure ADC level (arb)
    gen: gen_ch
    ro: ro_ch
    rolen: readout length (us)
    
    '''
    gench = np.array(['B', 'A'])
    roch = np.array(['D', 'C'])
    
    hw_cfg={"ro_ch" : ro,  # np.where(gench==ro)[0][0],   # ADC_D
        "res_ch"    : gen    # np.where(roch==gen)[0][0]    # DAC_B
       }
    
    readout_cfg={
        "readout_length":soccfg.us2cycles(rolen, gen_ch=gen), # [Clock ticks]
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
    phase_calibrated = phase_func(cal, phases, fpts)
    
    filename = '{}MHz to {}MHz pow {}adc_level'.format(fstart, fstop, pgain)+configs
    
    data = np.hstack((np.transpose([fpts]), np.transpose([amps]), np.transpose([phase_calibrated]), np.transpose([Is]), np.transpose([Qs])))
    np.savetxt(wdir+filename, data, delimiter=',', header='freq,Amp,Phase,I,Q\nMHz,V,Deg,V,V', comments='')
    
    # magindBm = (np.log(np.array(cal1[0])) - 5.637457678659047)/0.11333629477649314 - 2.3
        
    plt.figure(figsize=(20,15))
    plt.subplot(221,title="DUT", xlabel="Frequency (MHz)", ylabel="Amp. (adc level)")
    plt.plot(fpts, 20*np.log10(amps),'go-')
    max_freq=fpts[np.argmax(amps)]
    
    plt.subplot(222, title="phase", xlabel='Frequency', ylabel='Phase(degree)')
    plt.plot(fpts, phase_calibrated, '-')
    
    plt.subplot(212, projection='polar')
    #plt.axes(projection='polar')
    plt.plot(phase_calibrated*np.pi/180, amps, '.-')
    plt.grid(True)
    plt.show()
    
def resonator_spectroscopy(f_init, f_fin, f_step, temp, soccfg, soc, gen, ro, rolen, fstep, pgain, phcal, wdir):    
    
    flist = np.arange(f_init, f_fin, 0.8*f_step)
    
    for f in flist:
        f_start = f-f_step/2
        f_stop = f+f_step/2
        VNA(soccfg=soccfg, soc=soc, gen=gen, ro=ro, rolen=rolen, fstart=f_start, fstop=f_stop, fstep=fstep, pgain=pgain, 
            configs=temp+'_{}MHz.DAT'.format(f), cal=phcal, wdir=wdir) 
        
    list_temp = [x for x in os.listdir(wdir) if "DAT" and temp in x]
    freq_list = [x.split('_')[3].split('MHz')[0] for x in list_temp]
    
    list_temp_sort = zipsort(freq_list, list_temp)
    
    plt.figure(figsize=(10,5))
    plt.title(temp+' transmission')
    for dat in list_temp_sort[1]:
        spec = data_extract(Path=wdir, file=dat, data_type='rfsoc')
        plt.plot(spec[0], 20*np.log10(spec[1]), 'b-')
    plt.xlabel('Frequency(MHz)')
    plt.ylabel(r'Transmission $10*log_10(V_{adc}^2)$')  

def bfieldsweep(magnet, blo, bhi, bstep, ramptime, soccfg, soc, gen, ro, fstart, fstop, fstep, gain, temperature, inputpo, phcal, path):
    Bzlist1 = np.arange(0, bhi, bstep)
    Bzlist2 = np.arange(bhi, -blo, -bstep)
    Bzlist3 = np.arange(-blo, bhi, bstep)
    
    Bzlist = np.hstack((Bzlist1, Bzlist2, Bzlist3))
    
    #Bz.set_ramping(seg=1, rate=0.00005, upper_bound=1)
    
    for bz in Bzlist:
        
        magnet.set_field(bz)
        
        time.sleep(ramptime)
        
        VNA(soccfg=soccfg, soc=soc, gen=gen, ro=ro, fstart=fstart, fstop=fstop, fstep=fstep, pgain=gain, 
        configs='_{}mT_{}mK_input_{}dBm.DAT'.format(bz*1000, temperature, inputpo), cal=fitx, wdir=wdir_fieldsweep_RFSoC_new)

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

