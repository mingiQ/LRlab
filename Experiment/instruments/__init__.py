# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import Instrument, VisaInstrument, SocketInstrument, SerialInstrument
from Experiment.instruments.instrumentmanager import InstrumentManager


# try: from .InstrumentManagerWindow import InstrumentManagerWindow
# except: print("Could not load InstrumentManagerWindow")


from .TestInstrument import EchoInstrument,RandomInstrument
# from .TDS7104 import TekTDS7104
# from .RCA18 import MCRCA18
# from .multimeter import Keithley199

# try: from .relaybox.RFSwitch import RFSwitch
# except: print("Could not load heliummanifold")
# try: from .bkpowersupply import BKPowerSupply
# except: print("Could not load BKPowerSupply")
# try: from .bkpowersupply import BKPowerSupplynew
# except: print("Could not load BKPowerSupply")
# try: from .bkpowersupply import BKPowerSupply2
# except: print("Could not load BKPowerSupply")
try: from .mlBPF import MLBPF
except: print("Could not load MLBPF")
try: from .egg5210 import EGG5210
except: print("Could not load 5210")
try: from .anritsuMS2038 import MS2038
except: print("Could not load Anritsu")
try: from .temperature_Leiden import Thermo_leiden
except: print("Could not load leiden_thermo")
try: from .copperMT import CMT
except: print("Could not load Copper Mountain")
try: from .ami430magnet import AMI430
except: print("Could not AMI")
try: from .dmm2612 import Keithley_2612
except: print("Could not load Keithley 2612")
try: from .gigatronics1018 import GT1018
except: print("Could not load Gigatronics")
try: from .keithely2400 import Keithley_2400
except: print("Could not load K2400")
try: from .PNAX import PNA
except: print("Could not load PNAX")
try: from .lr_DAC import DAC_CF
except: print("Could not load DAC")
try: from .keithely237 import Keithley_237
except: print("Could not load Keithley 237")
try: from .rfgenerator import E4421B
except: print("Could not load E4421B")
try: from .spectrum_analyzer import E4440
except: print("Could not load E4440A")
try: from .aeroflex import Aeroflex
except: print("Could not load Aeroflex")
# try: from .KEPCOPowerSupply import KEPCOPowerSupply
# except: print("Could not load KEPCOPowerSupply")
# try: from .voltsource import SRS900
# except: print("Could not load SRS900")
# try: from .voltsource import YokogawaGS200
# except: print("Could not load YokogawaGS200")
# from .Alazar import Alazar, AlazarConfig, AlazarConstants
# try: from .Alazar import Alazar, AlazarConfig, AlazarConstants
# except: print("Could not load Alazar card")
# try: from .function_generator import BiasDriver,FilamentDriver,BNCAWG

# try: from .lockin import SR844
# except: print("Could not load SR844 driver")
# =============================================================================
# try: from .PressureGauge import PressureGauge
# except: print("Could not load PressureGauge driver")
# try: from .RGA100 import RGA100
# except: print('Could not load SRS RGA100 driver')
# try: from .AG850 import AG850
# except: print("Could not load AG850 Driver")
# try: from .Autonics import TM4
# except: print("Could not load Autonics TM4 Driver")
# try: from .Oerlikon import Center_Three
# except: print("Could not load Oerlikon Center Three Driver")
# try: from .TempScanner import HP34970A
# except: print("Could not load HP34970A Driver")
# try: from .PLC import FurnacePLC
# except: print("Could not load Furnace PLC")
# =============================================================================
