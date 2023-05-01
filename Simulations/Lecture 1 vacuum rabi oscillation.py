# In[necessary pkgs]

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from qutip import about, basis, destroy, mesolve, ptrace, qeye, tensor, wigner

# In[Jaynes-Cumming model]

'''
hbar = 1
'''

wc = 1.0 * 2 * np.pi  # cavity freq
wq = 1.0 * 2 * np.pi  # qubit freq
g = 0.05 * 2 * np.pi  # coupling rate
kappa = 0.005         # cavity dissipation
gamma = 0.05          # qubit dissipation
N = 15                # number of cavity fock states
n_th_a = 0.0          # avg number of thermal bath excitation
use_rwa = True

tlist = np.linspace(0, 25, 101)


# In[operators, Hamiltonian and initial state]

# initial state
psi0 = tensor(basis(N,0), basis(2,1))  # photon number zero / qubit |1>

# operators
a = tensor(destroy(N), qeye(2))   # cavity annihilation
sm = tensor(qeye(N), destroy(2))  # qubit annihilation

# Hamiltonian

if use_rwa:
    H = wc * a.dag() * a + wq * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())
    
else:
    wc * a.dag() * a + wq * sm.dag() * sm + g * (a.dag() + a) * (sm + sm.dag()) 
    

# In[dissipation]

c_ops = []

# cavity relaxation
rate = kappa * (1+ n_th_a)
if rate > 0.0:
    c_ops.append(np.sqrt(rate) * a)
    

# cavity excitation, if temp > 0
rate = kappa * n_th_a 
if rate > 0.0:
    c_ops.append(np.sqrt(rate) * a.dag())
    
# qubit relaxation
rate = gamma
if rate > 0.0:
    c_ops.append(np.sqrt(rate) * sm)
    
# In[solve]

output = mesolve(H, psi0, tlist, c_ops, [a.dag() * a, sm.dag() * sm, H])

# In[visualize]

n_c = output.expect[0]
n_a = output.expect[1]
E_system = output.expect[2]

fig, axes = plt.subplots(1, 1, figsize=(10, 6))

axes.plot(tlist, n_c, label="Cavity")
axes.plot(tlist, n_a, label='atom excited state')
#axes.plot(tlist, E_system, label='energy expectation')
axes.legend(loc=0)
axes.set_xlabel("time")
axes.set_ylabel("occupation probability")
axes.set_title("vacuum rabi oscillation")


# In[make helper function]

def Jaynes_Cumming(fc, fq, g, kappa, gamma, N, n_th_a, tlist, use_rwa):
    wc = fc * 2 * np.pi  # cavity freq
    wq = fq * 2 * np.pi  # qubit freq
    g = g * 2 * np.pi  # coupling rate
    #kappa = 0.005         # cavity dissipation
    #gamma = 0.05          # qubit dissipation
    #N = 15                # number of cavity fock states
    #n_th_a = 0.0          # avg number of thermal bath excitation
    #use_rwa = True
    
    #tlist = np.linspace(0, 25, 101)
    
    # initial state
    psi0 = tensor(basis(N,0), basis(2,1))  # photon number zero / qubit |1>
    
    # operators
    a = tensor(destroy(N), qeye(2))   # cavity annihilation
    sm = tensor(qeye(N), destroy(2))  # qubit annihilation
    
    # Hamiltonian
    
    if use_rwa:
        H = wc * a.dag() * a + wq * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())
        
    else:
        H = wc * a.dag() * a + wq * sm.dag() * sm + g * (a.dag() + a) * (sm + sm.dag()) 
        
    
    
    c_ops = []
    
    # cavity relaxation
    rate = kappa * (1+ n_th_a)
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * a)
        
    
    # cavity excitation, if temp > 0
    rate = kappa * n_th_a 
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * a.dag())
        
    # qubit relaxation
    rate = gamma
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * sm)
    
    output = mesolve(H, psi0, tlist, c_ops, [a.dag() * a, sm.dag() * sm, H])
    
    
    
    n_c = output.expect[0]
    n_a = output.expect[1]
    E_system = output.expect[2]
    
    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    
    axes.plot(tlist, n_c, label="Cavity")
    axes.plot(tlist, n_a, label='atom excited state')
    #axes.plot(tlist, E_system, label='energy expectation')
    axes.legend(loc=0)
    axes.set_xlabel("time")
    axes.set_ylabel("occupation probability")
    axes.set_title("vacuum rabi oscillation")

# In[]
tl = np.linspace(0,100,400)
Jaynes_Cumming(fc=1.4, fq=1, g=0.05, kappa=0.005, gamma=0.05, N=15, n_th_a=0, tlist=np.linspace(0,25,101) , use_rwa=True)

# In[cavity wigner function]
output = mesolve(H, psi0, tlist, c_ops, [])

'''
returns list of density matrices for the system

'''
print(output)


'''
look at the wigner functions at the point in time when atom is in its gnd state ( t=5,15, 25)

'''

# 1 find the indices of the density matrices for the time we are interested in 

t_idx = np.where([tlist == t for t in [0.0, 5.0, 10.0, 15.0, 20.0, 25.0]])[1]

print(tlist[t_idx])

# 1.5 get a list density matrices

rho_list = [output.states[i] for i in t_idx]

# loop over the lis t of density matrices

xvec = np.linspace(-3,3, 200)

fig, axes = plt.subplots(1, len(rho_list), sharex=True, figsize=(3*len(rho_list), 3))

for idx, rho in enumerate(rho_list):
    
    rho_cavity = ptrace(rho, 0) # partial trace, source: 0
    
    W = wigner(rho_cavity, xvec, xvec)
    
    axes[idx].contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-0.25, 0.25), cmap=plt.get_cmap("RdBu"))
    
    axes[idx].set_title(r"$t=%.1f$" % tlist[t_idx][idx])
    
    
    
    

# In[plot #2]

t_idx = np.where([tlist == t for t in [0.0, 5.0, 10, 15, 20, 25]])[1]

rho_list = [output.states[i] for i in t_idx]

fig_grid = (2, len(rho_list) * 2)

for idx, rho in enumerate(rho_list):
    
    rho_cavity = ptrace(rho, 0)
    W = wigner(rho_cavity, xvec, xvec)
    ax = plt.subplot2grid(fig_grid, (0, 2*idx), colspan=2)
    ax.contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-0.25, 0.25), cmap=plt.get_cmap("RdBu"))

# plot the cavity occupation probability in the ground state
ax = plt.subplot2grid(fig_grid, (1, 1), colspan=(fig_grid[1] - 2))
ax.plot(tlist, n_c, label="Cavity")
ax.plot(tlist, n_a, label="Atom excited state")
ax.legend()
ax.set_xlabel("Time")
ax.set_ylabel("Occupation probability");













