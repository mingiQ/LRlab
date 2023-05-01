# In[necessary pkgs]

import pyEPR as epr
import sys
import numpy as np
import matplotlib.pyplot as plt

# In[1]

print('Parsing unis:  1um =', 
      epr.parse_entry('1um', 'meters'), 'meters')

print(f"""For   L_J = 11 nH, the Josephson junction energy is
      E_J = {epr.calcs.Convert.Ej_from_Lj(11, 'nH', 'GHz'):.1f} GHz""")
      
h = 6.62607015e-34

phi0 = 2.067e-15

Lj = 11e-9

Ej = phi0**2/Lj/(4*np.pi**2)

print("Ej = {}J , {}GHz".format(Ej, Ej/h/1e9))

# In[load HFSS tutorial]

# Load Path temporarily just to find where the tutorial folder is
# return path_to_project
from pathlib import Path
path_to_project = Path(epr.__file__).parent.parent / '_example_files'
print(f'We will the example project located in\n {path_to_project}')

# In[load tutorial]

pinfo = epr.ProjectInfo(project_path = path_to_project, 
                         project_name = 'pyEPR_tutorial1',
                         design_name  = '1. single_transmon')