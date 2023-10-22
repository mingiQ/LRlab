# In[0]:pkgs

import sys
import time
import numpy as np

sys.path.append('Z:/Mingi Kim/python_codes/Experiment/test_insttype')
import VNA as vna 

# In[]

va = vna.MS2038()

# In[]

va.idn()

# In[]

va.start(12e9)

# In[]

va.sweep("single")

# In[]

va.points("3333")

# In[]
va.IFBW(2000)
va.points(2000)
time.sleep(1)

va.sweep("single")

# In[]

va.save_s2p('Z:/Mingi Kim/python_codes/Experiment/test_insttype/vna_test.s2p')