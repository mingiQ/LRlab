# In[0]:
import matplotlib.pyplot as plt
import numpy as np
from qutip import Qobj, about, basis, steadystate, coherent, displace, squeeze, coherent_dm, create, destroy, expect, fock, fock_dm, mesolve, qeye, sigmax, sigmay, sigmaz, tensor, wigner, thermal_dm, plot_wigner
  
# In[1]:
'''    
Qobj class --> representing quantum object (states &  operator)

'''

# In[creating and inspecting quantum objects]

q = Qobj([[1],[0]])

# In[inspection of qobj]

print('the dimension'+str(q.dims)+'\n'
      +'the shape'+str(q.shape)+'\n'
      +'data itself'+str(q.data)+'\n'
      +'dense matix representation'+str(q.full())+'\n'
      +str(q.isherm) +'\n'
      +str(q.type))

# In[using qobj instances for calculations]

sy = Qobj([[0,-1j],[1j,0]])
sz = Qobj([[1,0],[0,-1]])

print(sy)

print(sz)

H = 1.0 * sz + 0.1 * sy

print("Qubit Hamiltonian = \n")
print(H)

print(sy.dag())

print('trace = ' + str(H.tr()))

print('Eigenenergies = {} '.format(H.eigenenergies()))


# In[more features of Qobj]'

# In[states and operators]

N = 2 # num of states in the Hilbert space
n = 1 # the state that will be occupied

print(basis(N,n))

print(fock(N,n))

# Coherent state
coh = coherent(N=10, alpha=5.0)
alphacoh = destroy(10) * coh
print(coh, alphacoh)  # ?

# In[Density matrices]

dm = fock_dm(5, 2) # 5: hilbert space, 2=state that is occupied
prod = basis(5, 2) * basis(5, 2).dag()
print(dm, prod)

cdm = coherent_dm(N=3, alpha=1.0) # Coherent state
prod = coherent(N=3, alpha=1.0) * coherent(N=3, alpha=1.0).dag()
print(cdm, prod)

n = 1 # avg num of thermal photons
tdm = thermal_dm(5, n, 'analytic')
print(tdm)

'''
diag(n**i / (n+1)**i , i:range(N))
'''

numop = create(N=5) * destroy(N=5)

O = tdm * numop
#print(O.tr()) 

# In[operators]

print(destroy(N=3))

print(create(N=3))


# In[composite systems]

'''
Describing coupled quantum systems (qubit - cavity)

'''

sz1 = tensor(sigmaz(), qeye(2))
N=2
# flip excited first qubit and unaffect second qubit 
print(sz1)

psi1 = tensor(basis(N,1), basis(N,0))
psi2 = tensor(basis(N,0), basis(N,1))

print(sz1 * psi2 - psi2)

'''
coupled two-qubit Hamiltonian

'''
epsilon = [1.0, 1.0]
g = 0.1

sz1 = tensor(sigmaz(), qeye(2))
sz2 = tensor(qeye(2), sigmaz())
xx = tensor(sigmax(), sigmax())

H = epsilon[0] * sz1 + epsilon[1] * sz2 + g * xx

print(H)

'''
Jaynes-Cumming Hamiltonian for a qubit-cavity system

'''

wc = 1.0 # cavity freq
wq = 0.9 # qubit freq
g = 0.1 # coupling rate

# cavity mode operator
a = tensor(destroy(5), qeye(2))

# qubit operators
sz = tensor(qeye(5), sigmaz())   # sigma z
sm = tensor(qeye(5), destroy(2)) # sigma minus

H0 = wc * a.dag() * a - 0.5 * wq * sz + g * (a * sm.dag() + a.dag() * sm)

print(H0)

# In[Hamiltonian spectroscopy]
r = []
fd_list = np.linspace(0.5, 1.5, 100)
gamma = 0.05
kappa = 0.005

for fd in fd_list:                             # qubit pump frequencies.
    H = H0 
    
    num_a = a.dag() * a
    num_b = sm.dag() * sm
    
    H -= fd*num_a #+ fp*num_a
    H -= wc*num_b 
 
 
    c_ops = []
    c_ops.append(np.sqrt(gamma)*a)
    c_ops.append(np.sqrt(kappa)*sm)
    
    rho_ss = steadystate(H, c_ops)
    r.append(expect(num_a,rho_ss))
    
# In[plot]:
    
plt.plot(fd_list, r)

# In[unitary dynamics]
'''
unitary evolution og quantum system --> master-equation solve (mesolve)

data type -- odedata

e.g. Hamiltonian H = sigmax init state = |0>

'''
H = sigmax()

psi0 = basis(2,0)

tlist = np.linspace(0, 10, 100)

result = mesolve(H, psi0, tlist, [], [])

'''
def mesolve(H, rho0, tlist, c_ops=[], e_ops=[], args={}, options=None,
            progress_bar=None, _safe_mode=True):

collapse operators : dissipation
    
'''
# In[results]

print(result)

print(result.states[-1])

print(expect(sigmaz(), result.states[-1]))


plt.plot(tlist, expect(sigmaz(), result.states))

plt.plot(tlist, expect(H, result.states))

'''
if we only interested in expectation values,

'''
result_exp = mesolve(H, psi0, tlist, [], [sigmax(), sigmay(), sigmaz()])


fig, axes = plt.subplots(1,1)

axes.plot(tlist, result_exp.expect[2], label=r"$\left< \sigma_z \right> $")
axes.plot(tlist, result_exp.expect[1], label=r"$\left< \sigma_y \right> $")
axes.plot(tlist, result_exp.expect[0], label=r"$\left< \sigma_x \right> $")

axes.set_xlabel(r"$t$", fontsize=20)
axes.legend(loc=2)

# In[dissipative dynamics]

'''
QHO loses photons to its env with a relaxation rate k

'''

w = 1.0         # oscillator freq
kappa = 0.1     # relaxation rate
a = destroy(100) # oscillator annhilation op
rho0 = fock_dm(100,99)  # initial state, fock state with 5 photons
H = w * a.dag() * a  # QHO

# a list of collapse operators
c_ops = [np.sqrt(kappa) * a]

# In[solve]

tlist = np.linspace(0, 50, 100)

'''
request that the solve return the expectation value
of the photon number state operator a.dag() * a

'''

result = mesolve(H, rho0, tlist, c_ops, [a.dag() * a])

# In[plot]

fig, axes = plt.subplots(1,1)

axes.plot(tlist, result.expect[0])
axes.set_xlabel(r"$t$", fontsize=20)
axes.set_ylabel(r"Photon number")


# In[wigner functions]

import matplotlib as mpl
import matplotlib.transforms as transforms
from ipywidgets import interact

# In[helper functions]

def plot_wigner_psi_phi(psi, alpha_max=7.5):
    
    fig = plt.figure(figsize=(9,9))
    widths = [6,3]
    heights = [6,3]
    spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=widths, height_ratios=heights)
    
    x = np.linspace(-alpha_max, alpha_max, 200)
    wig = wigner(psi, x, x)
    psi_x = np.sum(wig, axis=0)
    psi_p = np.sum(wig, axis=1)
    
    ax = fig.add_subplot(spec[0,0])
    plot_wigner(psi, fig=fig, ax=ax, alpha_max=alpha_max)
    
    ax = fig.add_subplot(spec[0,1])
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(90)
    
    ax.plot(x, -psi_p, transform=rot+base)
    ax.set_xticks([])
    ax.set_ylim(-alpha_max, alpha_max)
    
    ax = fig.add_subplot(spec[1,0])
    ax.plot(x, psi_x)
    ax.set_yticks([])
    ax.set_xlim(-alpha_max, alpha_max)
    
    
    
# In[Harmonic oscillators]

N = 10
psi = fock(N,0)

plot_wigner_psi_phi(psi)  

# In[coherent state]

N = 25
alpha = 2+2j
psi = coherent(N,alpha/np.sqrt(2))

plot_wigner_psi_phi(psi, alpha_max=7)  

# In[cat state]

N = 30
psi = coherent(N,-3j) + coherent(N,3j)
plot_wigner_psi_phi(psi,alpha_max=7)

# In[interactive]

N = 30
def update(phi=0):
    psi = coherent(N, 3*np.exp(1j*phi)) + coherent(N, -3*np.exp(1j*phi))
    plot_wigner_psi_phi(psi)
    
interact(update, phi=(0,2*np.pi, 2*np.pi/24))

# In[]

N = 30
psi = displace(N,4/np.sqrt(2)) * basis(N, 0) 
plot_wigner_psi_phi(psi)

# In[]

N = 30
psi = displace(N,4/np.sqrt(2)) * squeeze(N,3) * basis(N, 0) 
plot_wigner_psi_phi(psi)

# In[hw5]


def W_new(x,p):
    return x*np.exp(-x**2)*np.exp(-p**2)

xvec=np.linspace(-10,10,100)

fig = plt.figure(figsize=(9,9))
widths = [6,3]
heights = [6,3]
spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=widths, height_ratios=heights)
    
    
W = np.zeros((100,100))

    
for x in range(len(xvec)):
    for p in range(len(xvec)):
        W[x][p] = W_new(xvec[x], xvec[p])

psi_x = np.sum(W, axis=0)
psi_p = np.sum(W, axis=1)

ax = fig.add_subplot(spec[0,0])
ax.contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-0.25, 0.25), cmap=plt.get_cmap("RdBu"))

ax = fig.add_subplot(spec[0,1])
base = plt.gca().transData
rot = transforms.Affine2D().rotate_deg(90)
    
ax.plot(xvec, -psi_p, transform=rot+base)

ax = fig.add_subplot(spec[1,0])
ax.plot(xvec, psi_x)
ax.set_yticks([])

# In[quadrature operators]
x_mean = 0
x_sqmean = 0

for x in range(len(xvec)):
    x_mean += xvec[x]*psi_x[x]
    x_sqmean += xvec[x]**2*psi_x[x]
    
Dx = x_sqmean - x_mean**2
print(Dx)

p_mean = 0
p_sqmean = 0

for p in range(len(xvec)):
    p_mean += xvec[p]*psi_p[p]
    p_sqmean += xvec[p]**2*psi_p[p]
    
Dp = p_sqmean - p_mean**2
print(Dp)
    

















