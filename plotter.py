import numpy as np
from CCD_pairing import *
from hamiltonian import *
import matplotlib.pyplot as plt
import sys

particles = 4
pairs = 2

g = 1
delta = 1

if g > 0:
    alpha = 0.8
else:
    alpha = 0.4

n = 20 # Basis
 
CI_basis = np.array(list(range(particles,n+pairs)))
CC_basis = CI_basis*2 - particles

length = CI_basis.shape[0]

CC1 = np.zeros(length)
CC1corr = np.zeros(length)
CI1 = np.zeros(length)
CI1corr = np.zeros(length)

for i in range(CI_basis.shape[0]):
    print('{} / {}'.format(i+1,CI_basis.shape[0]))
    b = CC_basis[i]
    E,E_ref = run(particles,b,delta,g,alpha=alpha)
    CC1[i] = E_ref+E
    CC1corr[i] = E
    H,_ = hamiltonian(pairs,CI_basis[i],delta,g)
    e,v = np.linalg.eigh(H)
    CI1[i] = e[0]
    CI1corr[i] = e[0]-H[0][0]


plt.rc('figure', titlesize=16)
plt.rc('axes', titlesize=14)
plt.rc('axes', labelsize=14)

x = CC_basis+particles

figsize = [7,5]

plt.figure(figsize=figsize)
plt.title('Ground state energy (g = {})'.format(1))
plt.plot(CC_basis+particles,CC1,label='CCD')
plt.plot(CC_basis+particles,CI1,label='FCI')
plt.xlabel('Number of spin orbitals')
plt.ylabel('Energy [Hartree]')
plt.legend()
plt.savefig('figures/4partE0_g.pdf')

plt.figure(figsize=figsize)
plt.title('Energy difference (g = {})'.format(1))
plt.plot(CC_basis+particles,np.abs(CC1-CI1),label='CCD')
plt.xlabel('Number of spin orbitals')
plt.ylabel('Energy [Hartree]')
plt.legend()
plt.savefig('figures/4partEdiff_g.pdf')

plt.figure(figsize = figsize)
plt.plot(CC_basis+particles,CC1corr,label='CCD')
plt.plot(CC_basis+particles,CI1corr,label='FCI')
plt.title('Correlation energy g=1')
plt.legend()
plt.xlabel('Number of spin orbitals')
plt.ylabel('Energy [Hartree]')

plt.show()
