# In[0] : necessary packages
import sys
import time
import numpy as np
from math import sqrt, ceil
sys.path.append('Z:/general/LRlabcode/LRlab')
import pyvisa
from Experiment.instruments.instrumenttypes import VisaInstrument

DEFAULT_TIME_STEP = 0  # in seconds
DEFAULT_NUM_POINTS = 1  # number of data points to collect for each measurement

### Keithley digital multimeter 2400   :  https://github.com/hinnefe2/keithley/blob/master/keithley.py
class Keithley_2400:
    
    def __init__(self, name='K2612', address='GPIB0::19::INSTR', enabled=True, timeout=0.1):
        VisaInstrument.__init__(self, name, address, enabled)
        self.query_sleep=0.05
        self.timeout = 160000
        self.term_char=''
        
    def get_idn(self):
        return(self.instrument.query("*IDN?"))


    #####################################################################################################
    # Internal methods: these are used internally but shouldn't be necessary for basic use of the class #
    #####################################################################################################

    # do setup stuff I don't really understand
    # adapted from http://pyvisa.sourceforge.net/pyvisa.html#a-more-complex-example
    
    
    
    def _initialize(self):
        self.instrument.write("*RST")
        self.instrument.write("*CLS")
        self.instrument.write("STATUS:MEASUREMENT:ENABLE 512")
        self.instrument.write("*SRE 1")
        self.instrument.write("ARM:COUNT 1")
        self.instrument.write("ARM:SOURCE BUS")
        self.instrument.write("TRACE:FEED SENSE1")
        self.instrument.write("SYSTEM:TIME:RESET:AUTO 0")
        
        # set various things to default values
        #self.setDelay()
        self.setNumPoints()
    
    
    # clear the saved data from previous measurement
    def _clearData(self):
        self.dataAll = []
        self.dataVolt = []
        self.dataCurr = []
        self.dataRes = []
        self.dataTime = []
        self.data = {}
    
    
    # start a measurement and wait for the 'measurement is done' signal from the Keithley
    def _startMeasurement(self):
        self._startNoWait()
        #self._catchSRQ()
    
    
    # start a measurement
    def _startNoWait(self):
        self.instrument.write("OUTPUT ON")
        self.instrument.write("TRACE:FEED:CONTROL NEXT")    # "trace:clear" cannot clear the buffer in this mode
        self.instrument.write("INIT")
        #self.instrument.write("*TRG")
        #self.instrument.trigger()
        
    
    def _traceclear(self):
        self.instrument.write("TRACE:FEED:CONTROL NEVer") 
        self.instrument.write("TRACE:CLEAR")
        
    # stop a measurement, turn output off
    def _stopMeasurement(self):
        self.instrument.write("OUTPUT OFF")
        self.instrument.write("TRACE:FEED:CONTROL NEVer") 
        self.instrument.write("TRACE:CLEAR")
        #self.instrument.query("STATUS:MEASUREMENT?")
        
    # pull data from the Keithley
    # always call this before _stopMeasurement() bc _stopMeasurement clears the keithley's buffer
    def _pullData(self, printresult=False):
        # returns (V, I, I/V, time, status) for each data point
        # (at least when measuring resistance)
        # when not measuring resistance, I/V column = 9.91e37
        #self.instrument.timeout = 1000
        
        #+1.000206E+00, +1.000000E-04, +1.000236E+04, +7.282600E+01, +4.813200E+04
        try: 
            datastring = self.instrument.query("TRACE:DATA?")
            data = datastring.split('\n')[0].split(',')
            #dataarray = np.array([float(x) for x in data]).reshape(len(data)//5, 5)
            dataarray = np.array([float(x) for x in data])
        
        except pyvisa.VisaIOError:
            dataarray = 'no data in buffer!'
        
               
        return dataarray

        
     # get what is being measured (VOLTage or CURRent or RESistance)
    def getMeasure(self):
        # keithley returns something like ' "VOLT:DC", "RES" ' or ' "CURR:DC" '
        return self.instrument.query("SENSE:FUNCTION?").split('\n')[0].strip('""')

    # get what is being sourced (VOLTage or CURRent)
    # returns ['VOLTAGE'|'CURRENT', value in volts|amps]
    def getSource(self):
        source = self.instrument.query("SOURCE:FUNCTION:MODE?")
        
        if source == "VOLT\n":
            source_level = self.instrument.query("SOURCE:VOLTAGE:LEVEL?")
            return ['volt', float(source_level.split('_')[0])]
        elif source == "CURR\n":
            source_level = self.instrument.query("SOURCE:CURRENT:LEVEL?")
            return ['current', float(source_level.split('_')[0])]
        
        
        
     
    ##############################################################
    # Configuration methods: use these to configure the Keithley #
    ##############################################################

    # set the number of data points to take
    def setNumPoints(self, numPts=DEFAULT_NUM_POINTS):
        self.instrument.write("TRIGGER:COUNT %d" % numPts)
        self.instrument.write("TRACE:POINTS %d" % numPts)

    # set the delay between data points (in sec)
    def setDelay(self, delay=DEFAULT_TIME_STEP):
        self.instrument.write("TRIGGER:DELAY %f" % delay)

    # set DC source, expects source to be either "voltage" or "current"
    def setSourceDC(self, source, value=0):
        
        if source.lower() == "voltage":
            self.instrument.write("SOURCE:FUNCTION:MODE VOLTAGE")
            self.instrument.write("SOURCE:VOLTAGE:MODE FIXED")
            self.instrument.write("SOURCE:VOLTAGE:RANGE " + str(value))
            self.instrument.write("SOURCE:VOLTAGE:LEVEL " + str(value))
        elif source.lower() == "current":
            self.instrument.write("SOURCE:FUNCTION:MODE CURRENT")
            self.instrument.write("SOURCE:CURRENT:MODE FIXED")
            self.instrument.write("SOURCE:CURRENT:RANGE " + str(value))
            self.instrument.write("SOURCE:CURRENT:LEVEL " + str(value))
            
        elif (self.getMeasure()=='RES'):
            self.instrument.write("SENSE:RESISTANCE:MODE MANUAL")

    # set sweep source, expects source to be either "voltage" or "current"
    def setSourceSweep(self, source, startValue, stopValue, sourceStep, timeStep=DEFAULT_TIME_STEP):
        numPts = ceil(abs((stopValue - startValue) / sourceStep)) + 1
        if self.getMeasure() == 'RES':
            self.instrument.write("SENSE:RESISTANCE:MODE MANUAL")
        if source.lower() == "voltage":
            self.instrument.write("SOURCE:FUNCTION:MODE VOLTAGE")
            self.instrument.write("SOURCE:VOLTAGE:MODE SWEEP")
            # self.write("SOURCE:VOLTAGE:RANGE " + str(stopValue))
            self.instrument.write("SOURCE:VOLTAGE:START " + str(startValue))
            self.instrument.write("SOURCE:VOLTAGE:STOP " + str(stopValue))
            self.instrument.write("SOURCE:VOLTAGE:STEP " + str(sourceStep))
            self.setNumPoints(numPts)
            self.setDelay(timeStep)
        elif source.lower() == "current":
            self.instrument.write("SOURCE:FUNCTION:MODE CURRENT")
            self.instrument.write("SOURCE:CURRENT:MODE SWEEP")
            self.instrument.write("SOURCE:CURRENT:RANGE " + str(stopValue))
            self.instrument.write("SOURCE:CURRENT:START " + str(startValue))
            self.instrument.write("SOURCE:CURRENT:STOP " + str(stopValue))
            self.instrument.write("SOURCE:CURRENT:STEP " + str(sourceStep))
            self.setNumPoints(numPts)
            self.setDelay(timeStep)
        else:
            print("Error: bad arguments")
        return numPts

    # set what is being measured (VOLTage or CURRent or RESistance)
    def setMeasure(self, measure):
        self.instrument.write("SENSE:FUNCTION:OFF 'CURR:DC', 'VOLT:DC', 'RES'")
        if measure.lower() == "voltage":
            self.instrument.write("SENSE:FUNCTION:ON 'VOLTAGE:DC'")
        elif measure.lower() == "current":
            self.instrument.write("SENSE:FUNCTION:ON 'CURRENT:DC'")
        elif measure.lower() == "resistance":
            self.instrument.write("SENSE:FUNCTION:ON 'CURRENT:DC'")
            self.instrument.write("SENSE:FUNCTION:ON 'RESISTANCE'")
        else:
            print("Expected one of [current, voltage, or resistance]")

    # set the upper limit for how much current / voltage will be sourced
    def setCompliance(self, source, limit):
        if source.lower() == 'voltage':
            self.instrument.write("SENS:VOLT:PROT " + str(limit))
        if source.lower() == 'current':
            self.instrument.write("SENS:CURR:PROT " + str(limit))
        else:
            print("Expected one of [current, voltage]")

    # set resistance measurements to 4-wire
    def setFourWire(self):
        if self.instrument.query("SENSE:FUNCTION?").split(",")[-1].strip('"') == 'RES':
            self.instrument.write("SYSTEM:RSENSE ON")
            print('Resistance measurement changed to 4-wire')
        else:
            print('Must be measuring resistance')

    # set resistance measurements to 2-wire
    def setTwoWire(self):
        if self.instrument.query("SENSE:FUNCTION?").split(",")[-1].strip('"') == 'RES':
            self.instrument.write("SYSTEM:RSENSE OFF")
            print('Resistance measurement changed to 2-wire')
        else:
            print('Must be measuring resistance')

    # set triggering to use TLINK connections (for fastest linking of two Keithleys)
    def setTLINK(self, inputTrigs, outputTrigs):
        self.instrument.write("TRIG:SOURCE TLINK")
        self.instrument.write("TRIG:INPUT {}".format(inputTrigs))
        self.instrument.write("TRIG:OUTPUT {}".format(outputTrigs))

    # set triggering to be immediate (default for single Keithley measurements)
    def setNoTLINK(self):
        self.instrument.write("TRIG:SOURCE IMMEDIATE")
        self.instrument.write("TRIG:INPUT NONE")
        self.instrument.write("TRIG:OUTPUT NONE")


    ########################################################
    # Operation methods: use these to operate the Keithley #
    ########################################################

    # turn the output on
    def outputOn(self):
        self.instrument.write("OUTPUT ON")

    # turn the output off
    def outputOff(self):
        self.instrument.write("OUTPUT OFF")
        
        
# =============================================================================
# # In[]
# sm = Keithley_2400()
# 
# 
# # In[]
# 
# a = sm._pullData()
# 
# print(a)
# 
# # In[]
# 
# #sm.setCompliance('voltage', 5)
# sm.setSourceDC('current', 0.020)
# 
# sm._startNoWait()
# 
# # In[]
# 
# sm._stopMeasurement()
# # In[]
# import pyvisa
# 
# RM =  pyvisa.ResourceManager()
# sm = RM.open_resource('GPIB0::19::INSTR') 
# 
# sm.write("trace:clear")
# 
# # In[]
# 
# #a = sm.query("SOURCE:FUNCTION:MODE?")
# 
# #a = sm.ask_for_values("SOURCE:CURRENT:LEVEL?")
# 
# #a = sm.query("sense:function?")
# 
# a = sm.query("TRACE:DATA?")
# print(a)
# 
# # In[]
# 
# class Keithley2400(GpibInstrument):
#     """A class to interface with the Keithley 2400 sourcemeter"""
# 
#     def __init__(self, GPIBaddr):
#         try:
#             # call the visa.GpibInstrument init method w/ appropriate argument
#             super(Keithley2400, self).__init__("GPIB::%d" % GPIBaddr)
#             self._initialize()
#             self._clearData()
#             self.saveCounter = 0
#         except VisaIOError:
#             print('VisaIOError - is the keithley turned on?')
# 
#     # read a single data point
#     def measurePoint(self):
#         self.write("TRACE:FEED:CONTROL NEXT")
#         self.write("INIT")
#         self.trigger()
#         return self._pullData()
# 	
# 
#     # perform a measurement w/ current parameters
#     def doMeasurement(self):
#         self._clearData()
#         self._startMeasurement()
#         self._pullData()
#         self._stopMeasurement()
# 
#     # ramp the output from rampStart to rampTarget
#     def rampOutput(self, rampStart, rampTarget, step, timeStep=50E-3):
#         if rampTarget < rampStart: step = -abs(step)
#         if rampTarget > rampStart: step = abs(step)
# 
#         source = self.getSource()[0]  # either 'voltage' or 'current'
#         sourceValue = rampStart
#         self.setSourceDC(source, sourceValue)
#         self.write("OUTPUT ON")
#         # hack-y : while (current level + step) is closer to target than current level
#         while abs((sourceValue + step) - rampTarget) <= abs(sourceValue - rampTarget):
#             sourceValue += step
#             self.setSourceDC(source, sourceValue)
#             time.sleep(timeStep)
#         return sourceValue
# 
#     # starting with the output off, turn the output on then ramp the output up/down to a specified level
#     def rampOutputOn(self, rampTarget, step, timeStep=50E-3):
#         rampStart = 0
#         sourceValue = self.rampOutput(rampStart, rampTarget, step, timeStep)
#         return sourceValue
# 
#     # starting with the output on, ramp the output to 0, then turn the output off
#     def rampOutputOff(self, rampStart, step, timeStep=50E-3):
#         rampTarget = 0
#         sourceValue = self.rampOutput(rampStart, rampTarget, step, timeStep)
#         self.outputOff()
#         return sourceValue
# 
#     # save the collected data to file
#     # mode 'a' appends to existing file, mode 'i' increments file counter ie test0001.txt, test0002,txt
#     def saveData(self, filePath=DEFAULT_SAVE_PATH, fileName="test.txt", mode='i'):
#         if filePath[-1] != '/': filePath += '/'
#         if fileName[-4:] != '.txt': fileName += '.txt'
# 
#         if mode == 'a':
#             saveFile = open(filePath + fileName, "a+")
#         elif mode == "i":
#             self.saveCounter = 0
#             while True:
#                 self.saveCounter += 1
#                 # print("checking file: " + filePath + fileName[:-4] + "{:04d}".format(self.saveCounter) + ".txt")
#                 # print("")
#                 if not os.path.exists(filePath + fileName[:-4] + "{:04d}".format(self.saveCounter) + ".txt"):
#                     break
#             saveFile = open(filePath + fileName[:-4] + "{:04d}".format(self.saveCounter) + ".txt", "a+")
#         else:
#             print("invalid mode")
#             return -1
# 
#         saveFile.write("\n")
#         # single line format
#         # saveFile.write(DEFAULT_ROW_FORMAT_HEADER.format("V (volts)","I (amps)","I/V (ohms)","t (s)","?"))
# 
#         # format for importing to Origin
#         saveFile.write(DEFAULT_ROW_FORMAT_HEADER.format("V", "I", "I/V", "t", "?"))
#         saveFile.write("\n")
#         saveFile.write(DEFAULT_ROW_FORMAT_HEADER.format("volts", "amps", "ohms", "s", "?"))
#         saveFile.write("\n")
#         for row in chunks(self.dataAll, 5):
#             saveFile.write(DEFAULT_ROW_FORMAT_DATA.format(*row))
#             saveFile.write("\n")
#         saveFile.close()
# 	return saveFile.name
# 
#     def printSummary(self):
#         print("Measuring: " + self.getMeasure())
#         print("Sourcing: " + str(self.getSource()))
#         print("")
#         print(DEFAULT_ROW_FORMAT_HEADER.format("V (volts)", "I (amps)", "I/V (ohms)", "t (s)", "?"))
#         for row in chunks(self.dataAll, 5):
#             print(DEFAULT_ROW_FORMAT_DATA.format(*row))
#         print("")
# 
# 
# 
# 
# =============================================================================
