# In[0] : necessary packages
import numpy as np
import pyvisa
import serial
import socket
import json
import time
import sys

## DAC LR 
class DAC_CF:
    '''
    !!! different from BPF, PySerial requires '\n' at the end of the command!!!!! (2023/5/9)
        /*    ------ communication commands ---------
     *     serial @ 115200, line end character: \n
     *    "?" or "*IDN?"             ->  devstamp         - equipment ID  
     *    "DD dac"                   -> "DD"              - initialize dac
     *    "WF dac val val"           -> "WF"              - set single dac value (float)
     *    "WD dac dval"              -> "WD"              - set single dac value (digital), always in dc mode
     *    "WH dac val val"           -> "WH"              - set double dac & dac+1 to val
     *    "RF dac"                   -> "RF val"          - read dac value (float)
     *    "RD dac"                   -> "RD dval"         - read dac value (digital)
     *    "RH dac"                   -> "RH val"          - read double-dac value
     *    
     *    "WS dac val"               -> "WS"              - write dc change rate [V/s], 0=no ramp
     *    "RS dac"                   -> "RS val"          - read  dc change rate [V/s], 0=no ramp
     *    
     *    "WA dac ac                 -> "WA ac"           - set amplitude of ac modulation in V
     *    "WT dac freq               -> "WT"              - set freq of ac modulation in [Hz]
     *    "WP dac phase              -> "WP"              - set phase of ac modulation in [deg]
     *    "WU dac duty               -> "WU"              - set pulse duty cycle (0-1)
     *    "WM dac mode               -> "WM"              - set modulation type 0-off, 1-sin, 2-square, 3-pulse
     *    "RA dac                    -> "RA float"        - read amplitude of ac modulation in Vrms
     *    "RT dac                    -> "RT float"        - read period of ac modulation in [Hz], 0-disabled
     *    "RP dac                    -> "RP float"        - read phase of ac modulation in [deg]
     *    "RU dac                    -> "RU float"        - read pulse duty cycle
     *    "RM dac                    -> "RM int"          - read modulation type 0-off, 1-sin, 2-square, 3-sin(fast)
     *    "RW dac                    -> "RW many"         - read modulation params "dc,ac,f,ph,du'
     *    
     *    "SR dac"                   -> "SR bits neg pos step dc ldac" - read dac info from ROM
     *    "SW dac neg pos bits"      -> "SW"              - write dac info to ROM
     *    "ER dac"                   -> "ER neg pos bits" - read dac info from EEPROM
     *    "SE dac"                   -> "SE"              - write data from ROM to EEPROM
     *    "ES dac"                   -> "ES"              - write data from EEPROM to ROM
     *    "EC #"                     -> "EC"              - TFT write echo level, 0 (off) 1 (+values) 2 (+command) 3 (+debug)
     *                                                          -1 supresses com responce for WF & WD
     *    
     *    "IC dac Vret Vmin Vmax Vthp Vthn t1 t2 t3 tramp repeats AIN" 
     *                               -> "IC Vres"         - set Ic sweep, returns V step size. AIN > 0(ADC) or 0 (comparator)
     *    "RC dac"                   -> "RC Vret Vmin Vmax Vthp Vthn t1 t2 t3 tramp repeats AIN Vres" 
     *    "IS dac mode"              -> "IS stat/stream"  - performs Ic sweeps, mode=0: <Ic>, std, min, max, repeats
     *                                                                               1: stream of Isw's
     *                                                                               2: strem (Vdac,Vadc) for each step between Vmin and Vmax
     *    "IT dac"                   -> "IT time1, time2" - returns max time in sec for n and 1 sweeps (for IS & IR commands)
     *    "SC AIN+, AIN- event"      -> "SC"              - enables  comparator CMP1
     *    "DC"                       -> "DC"              - disables comparator CMP1
     *    "DI"                       -> "DI"              - enable internal DAC0 (12bit, 3.3V)
     *    "DW val"                   -> "DW"              - set DAC0 val
     *    "PW val"                   -> "PW"              - set pin 23 0/1 LOW/HIGH
     *    
     *    "AR AIN"                   -> "AR val"          - reads value on ADC AIN
     *    "AS ave, cspeed, sspeed, res"  -> "AS rate"     - sets teensy ADC parameters, returns rate [us/read]
     *    
     *    "TS dac dval"              -> "TS"              - comm test: write xAAAAA to dac dval times
     *    "TX pin dval"              -> "TX"              - pin test: write HIGH/LOW to pin dval times
     *    "TT pin"                   -> "TT"              - toggle pin
     *    "ZZ"                       -> "ZZ txt"          - testing output 
     *    
     *    "IRQ"                                           - sets IRQready=false at a serial interrupt level
    '''
    def __init__(self, addr):       # COM4
        self.DAC = serial.Serial(port=str(addr), baudrate=115200, timeout=None)
        self.DAC.flushInput()
        self.DAC.flushOutput()
    def write_read(self, x):                    # CAUTION: all the command should end with '\n'
        self.DAC.write(bytes( x +'\n', 'utf-8'))
        time.sleep(0.05)
        #data = self.arduino.readline().decode('utf-8').rstrip()
        data = self.DAC.readline().rstrip()
        print(data)
    def idn(self):
        self.DAC.write(bytes('*IDN?\n', 'utf-8'))
        time.sleep(0.05)
        data = self.DAC.readline().rstrip()
        print(data)
    def voltsweep(self, dac_CF, volt):
        self.DAC.write(bytes("WF "+ str(dac_CF) +" "+ str(volt) +" "+ str(volt)+"\n", 'utf-8'))
        time.sleep(0.05)
        data = self.DAC.readline().rstrip()
        print(data)
    def quit_dac(self):
        self.DAC.close()