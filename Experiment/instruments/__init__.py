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
