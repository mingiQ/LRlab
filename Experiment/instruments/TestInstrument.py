# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 10:46:40 2023

@author: leogroup
"""

# In[0] : necessary packages
import sys
sys.path.append('Z:/general/LRlabcode/LRlab')
from Experiment.instruments.instrumenttypes import Instrument
import random


class EchoInstrument(Instrument):

    def  __init__(self,name,address='',enabled=True):
        Instrument.__init__(self,name,address='',enabled=True)
        self.s=''

    def write(self,s):
        self.s+=s
    
    def read(self):
        s=self.s
        self.s=''
        return s

    def echo(self,s):
        return s
        
class RandomInstrument(Instrument):
    def  __init__(self,name,address='',enabled=True):
        Instrument.__init__(self,name,address='',enabled=True)

    def read(self):
        return random.random()
        
    def random(self):
        return self.read()
