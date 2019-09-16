import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

CI_g1 = pd.read_csv('Data/FCITime_1.csv')
CI_gmin1 = pd.read_csv('Data/FCITime_min1.csv')
CC_g1 = pd.read_csv('Data/CCTime_1.csv')
CC_gmin1 = pd.read_csv('Data/CCTime_min1.csv')

# CC Alpha
CC_gmin1_a00 = pd.read_csv('Data/CC_g=-1_a=0.0.csv')
CC_gmin1_a02 = pd.read_csv('Data/CC_g=-1_a=0.2.csv')
CC_gmin1_a04 = pd.read_csv('Data/CC_g=-1_a=0.4.csv')
CC_gmin1_a06 = pd.read_csv('Data/CC_g=-1_a=0.6.csv')
CC_gmin1_a08 = pd.read_csv('Data/CC_g=-1_a=0.8.csv')
CC_gmin1_a10 = pd.read_csv('Data/CC_g=-1_a=1.0.csv')

#print(CI_g1)
#print(CC_g1)

CI_E_g1 = CI_g1['E']
CI_E_gmin1 = CI_gmin1['E']
CC_E_g1 = CC_g1['E']
CC_E_gmin1 = CC_gmin1['E']
#print(CC_gmin1_a00)
#print(CC_gmin1_a10)

x = CI_g1['number_of_basis_states']

#plt.figure()
#plt.title('Energy v number of basis states')
#plt.plot(x,CC_E_gmin1,label='CC')
#plt.plot(x,CI_E_gmin1,label='CI')
#plt.legend()

E00 = CC_gmin1_a00['E']
E02 = CC_gmin1_a02['E']
E04 = CC_gmin1_a04['E']
E06 = CC_gmin1_a06['E']
E08 = CC_gmin1_a08['E']
E10 = CC_gmin1_a10['E']

c00 = CC_gmin1_a00['conv']
c02 = CC_gmin1_a02['conv']
c04 = CC_gmin1_a04['conv']
c06 = CC_gmin1_a06['conv']
c08 = CC_gmin1_a08['conv']
c10 = CC_gmin1_a10['conv']

# ------- Energy Convergence ------- #
plt.figure()
plt.title('Energy convergence with alpha')
#plt.plot(x,E00,label=r'$\alpha=0.0$')
plt.plot(x,E02,label=r'$\alpha=0.2$')
plt.plot(x,E04,label=r'$\alpha=0.4$')
plt.plot(x,E06,label=r'$\alpha=0.6$')
plt.plot(x,E08,label=r'$\alpha=0.8$')
#plt.plot(x,E10,label=r'$\alpha=1.0$')
plt.plot(x,CI_E_gmin1,label='CI')
plt.legend()

# ------- Energy Difference ------- #
plt.figure()
plt.title(r'|$E_{CI} - E_{CC}$|')
plt.plot(x,np.abs(E02-CC_E_gmin1),label=r'$\alpha=0.2$')
plt.plot(x,np.abs(E04-CC_E_gmin1),label=r'$\alpha=0.4$')
plt.plot(x,np.abs(E06-CC_E_gmin1),label=r'$\alpha=0.6$')
plt.plot(x,np.abs(E08-CC_E_gmin1),label=r'$\alpha=0.8$')
plt.legend()

# ------- Iteration Convergence ------- #
plt.figure()
plt.title('convergence with alpha (Iterations)')
#plt.plot(x,c00,label=r'$\alpha=0.0$')
plt.plot(x,c02,label=r'$\alpha=0.2$')
plt.plot(x,c04,label=r'$\alpha=0.4$')
plt.plot(x,c06,label=r'$\alpha=0.6$')
plt.plot(x,c08,label=r'$\alpha=0.8$')
#plt.plot(x,c10,label=r'$\alpha=1.0$')
plt.legend()

plt.show()


