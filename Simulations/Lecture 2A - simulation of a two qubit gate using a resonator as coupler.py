# In[necessary pkgs]

import matplotlib.pyplot as plt
import numpy as np
from qutip import about, basis, concurrence, destroy, expect, fidelity, ket2dm, mesolve, phasegate, ptrace, qeye, sigmaz, sqrtiswap, tensor
from scipy.special import sici

# In[parameters]

N = 10

wc = 5.0 * 2 * np.pi
w1 = 3.0 * 2 * np.pi
w2 = 2.0 * 2 * np.pi

g1 = 0.01 * 2 * np.pi
g2 = 0.0125 * 2 * np.pi

tlist = np.linspace(0, 100, 500)

width = 0.5

# resonant SQRT iSWAP gate
T0_1 = 20
T_gate_1 = (1 * np.pi) / (4 * g1)

# resonant iSWAP gate
T0_2 = 60
T_gate_2 = (2 * np.pi) / (4 * g2)

# In[Operators, Hamiltonian and initial state]

'''
Kronecker prod : cavity (x) qubit1 (x) qubit2

'''

# cavity operators
a = tensor(destroy(N), qeye(2), qeye(2))
n = a.dag() * a

# qubit 1
sm1 = tensor(qeye(N), destroy(2), qeye(2))
sz1 = tensor(qeye(N), sigmaz(), qeye(2))
n1 = sm1.dag() * sm1

# qubit 2
sm2 = tensor(qeye(N), qeye(2), destroy(2))
sz2 = tensor(qeye(N), qeye(2), sigmaz())
n2 = sm2.dag() * sm2

# In[Hamiltonian]

Hc = a.dag() * a
H1 = -0.5 * sz1
H2 = -0.5 * sz2
Hc1 = g1 * (a.dag() * sm1 + a * sm1.dag())  # interaction cavity <-> qubit1
Hc2 = g2 * (a.dag() * sm2 + a * sm2.dag())  # interaction cavtiy <-> qubit2
H = wc * Hc + w1 * H1 + w2 * H2 + Hc1 + Hc2

print(H)

# In[initial state]
psi0 = tensor(basis(N,0), basis(2,1), basis(2,0))

# qubit 1 is in its excited state

# In[ideal two-qubit iSWAP gate]

def step_t(w1, w2, t0, width, t):
    '''
    step func that goes from w1 to w2 at time t0 as a fn of t
    
    '''
    return w1 + (w2-w1)*(t>t0)

fig, axes = plt.subplots(1,1, figsize=(8,2))
axes.plot(tlist, [step_t(0.5, 1.5, 50, 0.0, t) for t in tlist], 'k')
axes.set_ylim(0,2)
fig.tight_layout()

# In[modes]

def wc_t(t, args=None):
    return wc

def w1_t(t, args=None):
    return (w1 + step_t(0.0, wc-w1, T0_1, width, t)
                -step_t(0.0, wc-w1, T0_1 + T_gate_1, width, t))

def w2_t(t, args=None):
    return (w2 + step_t(0.0, wc-w2, T0_2, width, t)
                -step_t(0.0, wc-w2, T0_2 + T_gate_2, width, t))

H_t = [[Hc, wc_t], [H1, w1_t], [H2, w2_t], Hc1+Hc2]

# In[evolve the system]

res = mesolve(H_t, psi0, tlist, [], [])

# In[plot result]

fig, axes = plt.subplots(2,1, sharex=True, figsize=(12, 8))

axes[0].plot(tlist, np.array(list(map(wc_t, tlist)))/(2*np.pi), "r", linewidth=2, label='cavity')
axes[0].plot(tlist, np.array(list(map(w1_t, tlist)))/(2*np.pi), "b", linewidth=2, label='q1')
axes[0].plot(tlist, np.array(list(map(w2_t, tlist)))/(2*np.pi), "g", linewidth=2, label='q2')
axes[0].set_ylim(1,6)
axes[0].set_ylabel("energy(GHz)")
axes[0].legend()

axes[1].plot(tlist,  np.real(expect(n, res.states)), "r", label='cavity')
axes[1].plot(tlist,  np.real(expect(n1, res.states)), "b", label='q1')
axes[1].plot(tlist,  np.real(expect(n2, res.states)), "g", label='q2')
axes[1].set_xlabel("Time(ns)")
axes[1].set_ylabel("occupation probability")
axes[1].legend()

fig.tight_layout()

# In[inspect the final state]

rho_final = res.states[-1]

rho_qubits = ptrace(rho_final, [1,2])
print(rho_qubits)

# In[ideal iSWAP gate]
rho_qubits_ideal = ket2dm(tensor(phasegate(0), phasegate(-np.pi/2))*sqrtiswap() * tensor(basis(2,1), basis(2,0)))

print(rho_qubits_ideal) 

# In[fidelity and concurrence]

F = fidelity(rho_qubits, rho_qubits_ideal)
C = concurrence(rho_qubits) 

print(F, C)

# In[Finite pulse-rise time iSWAP gate]
def step_t(w1, w2, t0, width, t):
    return w1 + (w2-w1) / (1+np.exp(-(t-t0)/width))

fig, axes = plt.subplots(1, 1, figsize=(8, 2))
axes.plot(tlist, [step_t(0.5, 1.5, 50, width, t) for t in tlist], "k")
axes.set_ylim(0, 2)
fig.tight_layout()

# In[evolve the system]

res = mesolve(H_t, psi0, tlist, [], [])


# In[plot result]

fig, axes = plt.subplots(2,1, sharex=True, figsize=(12, 8))

axes[0].plot(tlist, np.array(list(map(wc_t, tlist)))/(2*np.pi), "r", linewidth=2, label='cavity')
axes[0].plot(tlist, np.array(list(map(w1_t, tlist)))/(2*np.pi), "b", linewidth=2, label='q1')
axes[0].plot(tlist, np.array(list(map(w2_t, tlist)))/(2*np.pi), "g", linewidth=2, label='q2')
axes[0].set_ylim(1,6)
axes[0].set_ylabel("energy(GHz)")
axes[0].legend()

axes[1].plot(tlist,  np.real(expect(n, res.states)), "r", label='cavity')
axes[1].plot(tlist,  np.real(expect(n1, res.states)), "b", label='q1')
axes[1].plot(tlist,  np.real(expect(n2, res.states)), "g", label='q2')
axes[1].set_xlabel("Time(ns)")
axes[1].set_ylabel("occupation probability")
axes[1].legend()

fig.tight_layout()

# In[Finite pulse-rise time & overshoot iSWAP gate]
def step_t(w1, w2, t0, width, t):
    return w1 + (w2-w1) * (0.5+sici((t-t0)/width)[0]/np.pi)

fig, axes = plt.subplots(1, 1, figsize=(8, 2))
axes.plot(tlist, [step_t(0.5, 1.5, 50, width, t) for t in tlist], "k")
axes.set_ylim(0, 2)
fig.tight_layout()

# In[evolve the system]

res = mesolve(H_t, psi0, tlist, [], [])


# In[plot result]

fig, axes = plt.subplots(2,1, sharex=True, figsize=(12, 8))

axes[0].plot(tlist, np.array(list(map(wc_t, tlist)))/(2*np.pi), "r", linewidth=2, label='cavity')
axes[0].plot(tlist, np.array(list(map(w1_t, tlist)))/(2*np.pi), "b", linewidth=2, label='q1')
axes[0].plot(tlist, np.array(list(map(w2_t, tlist)))/(2*np.pi), "g", linewidth=2, label='q2')
axes[0].set_ylim(1,6)
axes[0].set_ylabel("energy(GHz)")
axes[0].legend()

axes[1].plot(tlist,  np.real(expect(n, res.states)), "r", label='cavity')
axes[1].plot(tlist,  np.real(expect(n1, res.states)), "b", label='q1')
axes[1].plot(tlist,  np.real(expect(n2, res.states)), "g", label='q2')
axes[1].set_xlabel("Time(ns)")
axes[1].set_ylabel("occupation probability")
axes[1].legend()

fig.tight_layout()

# In[dissipative _ app collapse operators]

kappa = 0.00001  # high-Q
gamma1 = 0.005
gamma2 = 0.005

c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma1) * sm1, np.sqrt(gamma2) * sm2]

res = mesolve(H_t, psi0, tlist, c_ops, [] )

# In[plot result]

fig, axes = plt.subplots(2,1, sharex=True, figsize=(12, 8))

axes[0].plot(tlist, np.array(list(map(wc_t, tlist)))/(2*np.pi), "r", linewidth=2, label='cavity')
axes[0].plot(tlist, np.array(list(map(w1_t, tlist)))/(2*np.pi), "b", linewidth=2, label='q1')
axes[0].plot(tlist, np.array(list(map(w2_t, tlist)))/(2*np.pi), "g", linewidth=2, label='q2')
axes[0].set_ylim(1,6)
axes[0].set_ylabel("energy(GHz)")
axes[0].legend()

axes[1].plot(tlist,  np.real(expect(n, res.states)), "r", label='cavity')
axes[1].plot(tlist,  np.real(expect(n1, res.states)), "b", label='q1')
axes[1].plot(tlist,  np.real(expect(n2, res.states)), "g", label='q2')
axes[1].set_xlabel("Time(ns)")
axes[1].set_ylabel("occupation probability")
axes[1].legend()

fig.tight_layout()

# In[tunable coupler and fixed freq qubits]
width=0.15

def wc_t(t, arg=None):
    return (wc - step_t(0, wc-w1, T0_1, width, t)
                + step_t(0, wc-w1, T0_1+T_gate_1, width, t)
                - step_t(0, wc-w2, T0_2, width, t)
                + step_t(0, wc-w2, T0_2+T_gate_2, width, t))

def w1t(t):
    return w1

def w2t(t):
    return w2

H_t = [[Hc, wc_t], H1*w1 + H2*w2 + Hc1 + Hc2]

# In[evolve the system]

res = mesolve(H_t, psi0, tlist, c_ops, [])

# In[plot result]

fig, axes = plt.subplots(2,1, sharex=True, figsize=(12, 8))

axes[0].plot(tlist, np.array(list(map(wc_t, tlist)))/(2*np.pi), "r", linewidth=2, label='cavity')
axes[0].plot(tlist, np.array(list(map(w1t, tlist)))/(2*np.pi), "b", linewidth=2, label='q1')
axes[0].plot(tlist, np.array(list(map(w2t, tlist)))/(2*np.pi), "g", linewidth=2, label='q2')
axes[0].set_ylim(1,6)
axes[0].set_ylabel("energy(GHz)")
axes[0].legend()

axes[1].plot(tlist,  np.real(expect(n, res.states)), "r", label='cavity')
axes[1].plot(tlist,  np.real(expect(n1, res.states)), "b", label='q1')
axes[1].plot(tlist,  np.real(expect(n2, res.states)), "g", label='q2')
axes[1].set_xlabel("Time(ns)")
axes[1].set_ylabel("occupation probability")
axes[1].legend()

fig.tight_layout()

# In[tunable coupler and fixed freq qubits]
width=0.25

H_t = [[Hc, wc_t], [H1, w1_t], [H2, w2_t], Hc1+Hc2]

# In[evolve the system]

res = mesolve(H_t, psi0, tlist, c_ops, [])

# In[plot result]

fig, axes = plt.subplots(2,1, sharex=True, figsize=(12, 8))

axes[0].plot(tlist, np.array(list(map(wc_t, tlist)))/(2*np.pi), "r", linewidth=2, label='cavity')
axes[0].plot(tlist, np.array(list(map(w1_t, tlist)))/(2*np.pi), "b", linewidth=2, label='q1')
axes[0].plot(tlist, np.array(list(map(w2_t, tlist)))/(2*np.pi), "g", linewidth=2, label='q2')
axes[0].set_ylim(1,6)
axes[0].set_ylabel("energy(GHz)")
axes[0].legend()

axes[1].plot(tlist,  np.real(expect(n, res.states)), "r", label='cavity')
axes[1].plot(tlist,  np.real(expect(n1, res.states)), "b", label='q1')
axes[1].plot(tlist,  np.real(expect(n2, res.states)), "g", label='q2')
axes[1].set_xlabel("Time(ns)")
axes[1].set_ylabel("occupation probability")
axes[1].legend()

fig.tight_layout()