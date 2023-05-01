# In[0]: 
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from qutip import Options, about, basis, destroy, expect, mesolve, ptrace, qeye, sigmax, steadystate, tensor, wigner


# In[model]

'''
H = w0*a.dag()*a + 1/2*wq*sz + g*sx*(a.dag()+a)

kappa ~ cavity linewidth
Gamma ~ qubit excitation rate

dr/dt = -1j[H, r] + Gamma(sm.dag()*r*sm - 1/2*sm*sm.dag()*r - 1/2*r*sm*sm.dag())

kappa*(1+n_th)(a*r*a.dag() - 1.2*a.dag()*a*r - 1/2*r*a.dag()*a)
kappa*(n_th)(a.dag()*r*a - 1.2*a*a.dag()*r - 1/2*r*a*a.dag())

'''

w0 = 1.0 * 2 * np.pi  # cavity freq
wq = 1.0 * 2 * np.pi  # qubit freq
g = 0.05 * 2 * np.pi  # coupling rate

kappa = 0.04    # cavity linewidth
gamma = 0       # qubit dissipation
Gamma = 0.35    # atom pump rate

N = 50  # number of cavity fock states
n_th_a  = 0.0  # avg number pf thermal bath excitation

tlist = np.linspace(0, 150, 101)

# initial state
psi0 = tensor(basis(N,0), basis(2,0))

# operators
a = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
sx = tensor(qeye(N), sigmax())

# Hamiltoninan
H = w0*a.dag()*a + wq*sm.dag()*sm + g*(a.dag()+a)*sx

print(H)


# In[dissipation channel]

c_ops = []

rate = kappa * (1 + n_th_a)
if rate > 0:
    c_ops.append(np.sqrt(rate)*a)

rate = kappa * n_th_a
if rate > 0:
    c_ops.append(np.sqrt(rate)*a.dag())
    
rate = gamma
if rate > 0:
    c_ops.append(np.sqrt(rate)*sm)
    
rate=Gamma
if rate > 0:
    c_ops.append(np.sqrt(rate)*sm.dag()) 
    
# In[evolve the system]

opt = Options(nsteps=2000)  # allow extra time-steps
output = mesolve(H, psi0, tlist, c_ops, [a.dag()*a, sm.dag()*sm], options=opt)

# In[visualize the results]
n_c = output.expect[0]
n_q = output.expect[1]

fig, axes = plt.subplots(1,1, figsize=(8,6))

axes.plot(tlist, n_c, label='Cavity')
axes.plot(tlist, n_q, label='Atom excited state')
axes.legend(loc=0)
axes.set_xlabel("Time")
axes.set_ylabel("Occupation prob.")   

# In[steady state -- cavity fock-state dist. and wigner function]

rho_ss = steadystate(H, c_ops)

# In[plot]

fig, axes = plt.subplots(1, 2 , figsize=(12,6)) 

xvec = np.linspace(-5, 5, 200)

rho_cavity = ptrace(rho_ss, 0)
W = wigner(rho_cavity, xvec, xvec)
wlim = abs(W).max()

axes[1].contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-wlim, wlim), cmap=plt.get_cmap("RdBu"))
axes[1].set_xlabel(r"Im[$\alpha$]")
axes[1].set_ylabel(r"Re[$\alpha$]")

axes[0].bar(np.arange(0,N), np.real(rho_cavity.diag()), color='blue', alpha=0.6)
axes[0].set_ylim(0,1)
axes[0].set_xlim(0,N)
axes[0].set_xlabel("Fock number")
axes[0].set_ylabel("Occupation probability")















