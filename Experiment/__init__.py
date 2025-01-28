try:
    from .datamanagement import SlabFile, h5File, AttrDict, generate_file_path
except ModuleNotFoundError:
    print("one of the dependency is missing... ")

from .dataanalysis import *
from .matplotlib_text_wrapper import *
#from .experiment import Experiment
from . import kfit
