# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 10:41:58 2012

@author: Dave
"""

"""
To do this test you have to first start the nameserver 
and an InstrumentManager server

#To start a nameserver (if one isn't already present)
    python -m Pyro4.naming -n hostname
    
#To start InstrumentManager Server
    im = InstrumentManager(r'c:\_Lib\python\slab\instruments\instrument.cfg')
    im.serve_instruments()
    
"""
import sys
import time
sys.path.append('Z:/general/LRlabcode/LRlab')
import Experiment
from Experiment.instruments import InstrumentManager


im=InstrumentManager()


# In[]

print(list(im.keys()))
#print(im['echo'].echo('This is a test'))
#print(im['random'].random())
#print im['FRIDGE'].get_status()

# In[]

class C:
    def __init__(self, name, age):
        self.name = name
        self.age = age

    def m(self, x):
        print(f"{self.name} called with param '{x}'")
        return

ci = C("Joe", 10)

# In[]
print(C)
print(ci)
print(C.m)
print(ci.m)
print(getattr(ci,'m'))
getattr(ci,'m')('arg')

# In[]

import sys
class_module = sys.modules[Experiment.instruments.__class__.__module__]