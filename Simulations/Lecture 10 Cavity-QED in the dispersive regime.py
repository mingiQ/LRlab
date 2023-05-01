# In[0]: 
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from qutip import Options, about, basis, Bloch, coherent, destroy, expect, ket2dm, mesolve, ptrace, qeye, sigmax, sigmaz, steadystate, thermal_dm, tensor, correlation_2op_2t, spectrum_correlation_fft, wigner


# In[model]
'''
Dispersive regime occurs when the resonator and qubit is far off-resonance.
Delta >> g (Delta = wc-wq)

Dispersive regime 
H_int = g*(a.dag() + a)*sigmax() --> g**2/Delta * (num_cavity+1/2)*sigmaz()


'''

N = 20

wc = 2.0 * 2 * np.pi    # cavity mode
wq = 3.0 * 2 * np.pi    # qubit mode
chi = 0.025 * 2 * np.pi # parameter in the dispersive Hamiltonian

delta = abs(wc-wq)       # detuning
g = np.sqrt(delta * chi) # couping strength that is consistent with chi

print("dispersive regime: Delta = {} ,g = {}".format(delta/(2*np.pi) , g/(2*np.pi)))

# cavity operators
a = tensor(destroy(N), qeye(2))
nc = a.dag() * a
xc = a + a.dag()

# qubit operators
sm = tensor(qeye(N), destroy(2))
sz = tensor(qeye(N), sigmaz())
sx = tensor(qeye(N), sigmax())
nq = sm.dag() * sm
xq = sm + sm.dag()

Id = tensor(qeye(N), qeye(2))

# dispersive hamiltonian
H = wc*(a.dag()*a+Id/2) + (wq/2)*sz + chi * (a.dag()*a+Id/2) * sz

print(H)

# In[different initial states]

#psi0 = tensor(coherent(N, np.sqrt(6)), (basis(2,0)+basis(2,1)).unit())

#psi0 = tensor(coherent(N, np.sqrt(6)), basis(2,0))

#psi0 = tensor(coherent(N, np.sqrt(6)), basis(2,1))

psi0 = tensor(thermal_dm(N, 7), ket2dm(basis(2,0)+basis(2,1))).unit()


#psi0 = tensor(coherent(N, np.sqrt(4)), (basis(2, 0) + basis(2, 1)).unit() )



# In[evolve]

tlist = np.linspace(0, 250, 1000)
res = mesolve(H, psi0, tlist, [], [], options=Options(nsteps=5000))

# In[excitation numbers]
nc_list = expect(nc, res.states)
nq_list = expect(nq, res.states)

fig, ax = plt.subplots(1,1, sharex=True, figsize=(12, 4))

ax.plot(tlist, nc_list, "r", label="cavity")
ax.plot(tlist, nq_list, "b--", label="qubit")
ax.set_ylabel("n", fontsize=16)
ax.set_xlabel("Time (ns)", fontsize=16)
#ax.set_ylim(0, 6)
ax.legend()

fig.tight_layout()

# In[quadrature of the resonator & qubit]

xc_list = expect(xc, res.states)
xq_list = expect(xq, res.states)

fig, ax = plt.subplots(1,1, sharex=True, figsize=(12, 4))

ax.plot(tlist, xc_list, "r", label="cavity")
ax.plot(tlist, xq_list, "b--", label="qubit")
ax.set_ylabel("x", fontsize=16)
ax.set_xlabel("Time (ns)", fontsize=16)
#ax.set_ylim(0, 6)
ax.legend()

fig.tight_layout()

# In[correlation]

tlist = np.linspace(0, 100, 1000)
corr_vec = correlation_2op_2t(H, psi0, None, tlist, [], a.dag(), a)
corr_qubit = correlation_2op_2t(H, psi0, None, tlist, [], sx, sx)

# In[plot]

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12, 4))

ax.plot(tlist, np.real(corr_vec), "r", linewidth=2, label="resonator")
ax.plot(tlist, np.real(corr_qubit), "b--", linewidth=2, label="qubit")
ax.set_ylabel("correlation", fontsize=16)
ax.set_xlabel("Time (ns)", fontsize=16)
ax.legend()
ax.set_xlim(0, 50)
fig.tight_layout()

# In[spectral correlation function]

w, S = spectrum_correlation_fft(tlist, corr_vec)

# In[plot]

fig, ax = plt.subplots(figsize=(9,3))
ax.plot(w/2/np.pi, abs(S))
ax.set_xlabel(r"$\omega$")
ax.set_xlim(wc / (2 * np.pi) - 0.5, wc / (2 * np.pi) + 0.5)

# In[spectral correlation function]

w, S = spectrum_correlation_fft(tlist, corr_qubit)

# In[plot]

fig, ax = plt.subplots(figsize=(9,3))
ax.plot(w/2/np.pi, abs(S))
ax.set_xlabel(r"$\omega$")
#ax.set_xlim(wc / (2 * np.pi) - 0.5, wc / (2 * np.pi) + 0.5)

# In[]

fig, ax = plt.subplots(figsize=(9, 3))
ax.plot((w - wq - chi) / (2 * chi), abs(S), '.-')
ax.set_xlabel(r"$(\omega - \omega_q - \chi)/2\chi$", fontsize=18)
ax.set_xlim(-0.5, N);

# In[energy level]

def compute(w1list, wq, chi, N):

    # cavity operators
    a = tensor(destroy(N), qeye(2))
    
    # qubit operators
    sz = tensor(qeye(N), sigmaz())
    sx = tensor(qeye(N), sigmax())


    Id = tensor(qeye(N), qeye(2))
    
    # dispersive hamiltonian
    #H_kin = wc*(a.dag()*a+Id/2) + (wq/2)*sz 
    #H_int = g**2/Delta * (a.dag()*a+Id/2) * sz

    idx = 0
    evals_mat = np.zeros((len(w1list), 2*N))
    for w1 in w1list:

        # evaluate the Hamiltonian

        H_kin = w1*(a.dag()*a+Id/2) + (wq/2)*sz 
        H_int = chi * (a.dag() + a) * sx
        H = H_kin + H_int

        # find the energy eigenvalues of the composite system
        evals, ekets = H.eigenstates()

        evals_mat[idx, :] = np.real(evals)

        idx += 1

    return evals_mat

# In[eval]
w1 = 1.0 * 2 * np.pi  # atom 1 frequency: sweep this one
wq = 0.9 * 2 * np.pi  # atom 2 frequency
chi = 0.01 * 2 * np.pi  # cavity-qubit coupling strength

#w1list = np.linspace(0.75, 1.25, 50) * 2 * np.pi  # cavuty frequency range
w1list = np.linspace(0.75, 1.25, 100) * 2 * np.pi  # cavuty frequency range
evals_mat = compute(w1list, wq, chi, N=40)

# In[]
fig, ax = plt.subplots(figsize=(12, 12))

for n in range(40)[0:3]:
    ax.plot(w1list / (2*np.pi), (evals_mat[:, n] - 0*evals_mat[:, 0])/(2*np.pi), "-")

ax.set_xlabel("Energy splitting of cavity")
ax.set_ylabel("Eigenenergies")
ax.set_title("Energy spectrum of cavity coupled qubits")
ax.set_xlim(0.87, 1)

# In[energy level]

def compute_qubit(wqlist, w1, chi, N):

    # cavity operators
    a = tensor(destroy(N), qeye(2))
    
    # qubit operators
    sz = tensor(qeye(N), sigmaz())


    Id = tensor(qeye(N), qeye(2))
    
    # dispersive hamiltonian
    #H_kin = wc*(a.dag()*a+Id/2) + (wq/2)*sz 
    #H_int = g**2/Delta * (a.dag()*a+Id/2) * sz

    idx = 0
    evals_mat = np.zeros((len(w1list), 2*N))
    for wq in wqlist:

        H_kin = w1*(a.dag()*a+Id/2) + (wq/2)*sz 
        H_int = chi * (a.dag() + a) * sx
        H = H_kin + H_int

        # find the energy eigenvalues of the composite system
        evals, ekets = H.eigenstates()

        evals_mat[idx, :] = np.real(evals)

        idx += 1

    return evals_mat

# In[eval]
w1 = 1.0 * 2 * np.pi  # atom 1 frequency: sweep this one
wq = 0.9 * 2 * np.pi  # atom 2 frequency
chi = 0.04 * 2 * np.pi  # cavity-qubit coupling strength

w1list = np.linspace(0.7, 1.1, 50) * 2 * np.pi  # qubit frequency range

evals_mat = compute_qubit(w1list, w1, chi, N=20)

# In[]
fig, ax = plt.subplots(figsize=(12, 6))

for n in [0,1,2,3,4,5]:
    ax.plot(w1list / (2*np.pi), (evals_mat[:, n] - evals_mat[:, 0])/(2*np.pi), ".-")

ax.set_xlabel("Energy splitting of qubit")
ax.set_ylabel("Eigenenergies")
ax.set_title("Energy spectrum of cavity coupled qubits");

# In[energy level splitting]

def levelsplit(glist, wc, wq, N):

    # cavity operators
    a = tensor(destroy(N), qeye(2))
    
    # qubit operators
    sz = tensor(qeye(N), sigmaz())
    sx = tensor(qeye(N), sigmax())


    Id = tensor(qeye(N), qeye(2))
    
    # dispersive hamiltonian
    #H_kin = wc*(a.dag()*a+Id/2) + (wq/2)*sz 
    #H_int = g**2/Delta * (a.dag()*a+Id/2) * sz

    idx = 0
    evals_mat = np.zeros((len(w1list), 2*N))
    for g in glist:

        # evaluate the Hamiltonian

        H_kin = wc*(a.dag()*a+Id/2) + (wq/2)*sz 
        H_int = g * (a.dag() + a) * sx
        H = H_kin + H_int

        # find the energy eigenvalues of the composite system
        evals, ekets = H.eigenstates()

        evals_mat[idx, :] = np.real(evals)

        idx += 1

    return evals_mat
# In[eval]
w1 = 1.0 * 2 * np.pi  # atom 1 frequency: sweep this one
wq = 0.99 * 2 * np.pi  # atom 2 frequency
chi = 0.01 * 2 * np.pi  # cavity-qubit coupling strength

#w1list = np.linspace(0.75, 1.25, 50) * 2 * np.pi  # cavuty frequency range
glist = np.linspace(0, 5, 100) * 2 * np.pi  # cavuty frequency range
evals_mat = levelsplit(glist, w1, wq, N=40)

# In[]
fig, ax = plt.subplots(figsize=(10, 7))

for n in range(40)[2:10]:
    ax.plot(glist / (2*np.pi), (evals_mat[:, n] - 0*evals_mat[:, 0])/(2*np.pi), "-")

ax.set_xlabel("coupling rate")
ax.set_ylabel("Eigenenergies")
ax.set_title("Energy spectrum of cavity coupled qubits")
#ax.set_xlim(0.87, 1)

