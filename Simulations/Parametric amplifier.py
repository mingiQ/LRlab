# In[0]: 
import time

import matplotlib.pyplot as plt
import numpy as np
from qutip import (Odeoptions, about, basis, destroy, expect, mesolve, brmesolve,
                   propagator, propagator_steadystate, sigmax, sigmaz, plot_expectation_values)

# In[1]  Single mode - single port

N = 2

wa = 2 * 2*np.pi
gaa = 1 * 2*np.pi
Oa = 2 * wa

a = destroy(N)

H0 = wa * a.dag() * a 

# initial state
psi0 = basis(N, N - 1)

# times for simulation
times = np.linspace(0, 10, 100)

# solve using brmesolve
result_const = brmesolve(H0, psi0, times, e_ops=[a.dag() * a])

# In[expectation val]
plot_expectation_values(result_const, ylabels=["<n>"]);

# In[add time dependence]
def H1_coeff(t, arg):
    return np.sin(Oa*t)

Ht = [H0, [-gaa*a*a, H1_coeff],[-gaa*a.dag()*a.dag(), H1_coeff]]
result_me = mesolve(Ht, psi0, times, e_ops=[a.dag()*a])
plot_expectation_values(result_me, ylabels=["<n>"])

# In[]

# In[]

def hamiltonian_t(t, args):
    """evaluate the hamiltonian at time t."""
    H0 = args[0]
    H1 = args[1]

    return H0 + t * H1

def qubit_integrate(delta, eps0, A, gamma1, gamma2, psi0, tlist):

    # Hamiltonian
    sx = sigmax()
    sz = sigmaz()
    sm = destroy(2)

    H0 = -delta / 2.0 * sx - eps0 / 2.0 * sz
    H1 = -A / 2.0 * sz

    # collapse operators
    c_op_list = []

    n_th = 0.0  # zero temperature

    # relaxation
    rate = gamma1 * (1 + n_th)
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * sm)

    # excitation
    rate = gamma1 * n_th
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * sm.dag())

    # dephasing
    rate = gamma2
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * sz)

    # evolve and calculate expectation values

    # method 1: function callback which returns the time-depdent qobj
    # H_args = (H0, H1)
    # output = mesolve(hamiltonian_t,
    #                  psi0, tlist, c_op_list, [sm.dag() * sm], H_args)

    # method 2: a function callback that returns the coefficient for a qobj
    # H = [H0, [H1, lambda x,y: x]]
    # output = mesolve(H, psi0, tlist, c_op_list, [sm.dag() * sm], {})

    # method 3: a string that defines the coefficient. The solver generates
    # and compiles C code using cython. This method is usually the fastest
    # for large systems or long time evolutions, but there is fixed-time
    # overhead that makes it inefficient for small and short-time evolutions.
    H = [H0, [H1, "t"]]
    output = mesolve(H, psi0, tlist, c_op_list, [sm.dag() * sm], {})

    return output.expect[0]

# In[2]

# set up the calculation
#
delta = 0.5 * 2 * np.pi  # qubit sigma_x coefficient
eps0 = 0.0 * 2 * np.pi  # qubit sigma_z coefficient
A = 2.0 * 2 * np.pi  # sweep rate
gamma1 = 0.0  # relaxation rate
gamma2 = 0.0  # dephasing  rate
psi0 = basis(2, 0)  # initial state

tlist = np.linspace(-20.0, 20.0, 5000)
start_time = time.time()
p_ex = qubit_integrate(delta, eps0, A, gamma1, gamma2, psi0, tlist)
print("time elapsed = " + str(time.time() - start_time))

# In[3]

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(tlist, np.real(p_ex), "b", tlist, np.real(1 - p_ex), "r")
ax.plot(tlist, 1 - np.exp(-np.pi * delta**2 / (2 * A)) *
        np.ones(tlist.shape[0]), "k")
ax.set_xlabel("Time")
ax.set_ylabel("Occupation probability")
ax.set_title("Landau-Zener transition")
ax.legend(("Excited state", "Ground state", "Landau-Zener formula"), loc=0);