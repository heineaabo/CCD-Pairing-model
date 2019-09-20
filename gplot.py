import numpy as np
from matplotlib import pyplot as plt
import matplotlib.font_manager
import pandas as pd
import sys

save = False
if sys.argv[1] == 'true':
    save = True

##################################
########     g = 1.0      ########
##################################

#plt.rc('text',usetex=True)
font = {'family':'serif','size':14}
plt.rc('font',**font)

CI_g1 = pd.read_csv('Data/two_pairs/CI_g=1.csv')
CC_g1 = pd.read_csv('Data/CCTime_1.csv')

# CC Alpha
CC_g1_a02 = pd.read_csv('Data/two_pairs/CC_g=1_a=0.2.csv')
CC_g1_a04 = pd.read_csv('Data/two_pairs/CC_g=1_a=0.4.csv')
CC_g1_a06 = pd.read_csv('Data/two_pairs/CC_g=1_a=0.6.csv')
CC_g1_a08 = pd.read_csv('Data/two_pairs/CC_g=1_a=0.8.csv')
CC_g1_a10 = pd.read_csv('Data/two_pairs/CC_g=1_a=1.csv')

#print(CI_g1)
#print(CC_g1)

CI_E_g1 = CI_g1['E']
CI_t_g1 = CI_g1['Time']

x = CC_g1_a02['number_of_basis_states']

E02 = CC_g1_a02['E']
E04 = CC_g1_a04['E']
E06 = CC_g1_a06['E']
E08 = CC_g1_a08['E']
E10 = CC_g1_a10['E']

t02 = CC_g1_a02['Time']
t04 = CC_g1_a04['Time']
t06 = CC_g1_a06['Time']
t08 = CC_g1_a08['Time']
t10 = CC_g1_a10['Time']

c02 = CC_g1_a02['conv']
c04 = CC_g1_a04['conv']
c06 = CC_g1_a06['conv']
c08 = CC_g1_a08['conv']
c10 = CC_g1_a10['conv']

# ------- Energy Convergence ------- #
plt.figure(figsize=[7,5])
plt.title('Ground state energy',fontsize=20)
plt.xlabel('Number of basis states',fontsize=16)
plt.ylabel('Energy [Hartrees]',fontsize=16)
#plt.plot(x,E02,label=r'$\alpha=0.2$')
#plt.plot(x,E04,label=r'$\alpha=0.4$')
#plt.plot(x,E06,label=r'$\alpha=0.6$')
plt.plot(x,E08,label='CCD')
#plt.plot(x,E10,label=r'$\alpha=1.0$')
plt.plot(x,CI_E_g1,label='CI')
plt.legend()
if save:
    plt.savefig('figures/energy_g1_p2.pdf')

# ------- Energy Difference ------- #
plt.figure(figsize=[7,5])
plt.title(r'|$E_{CI} - E_{CC}$|',fontsize=20)
plt.xlabel('Number of basis states',fontsize=16)
plt.ylabel(r'$\Delta E$',fontsize=16)
#plt.plot(x,np.abs(E02-CI_E_g1),label=r'$\alpha=0.2$')
#plt.plot(x,np.abs(E04-CI_E_g1),label=r'$\alpha=0.4$')
#plt.plot(x,np.abs(E06-CI_E_g1),label=r'$\alpha=0.6$')
plt.plot(x,np.abs(E08-CI_E_g1))
#plt.plot(x,np.abs(E10-CI_E_g1),label=r'$\alpha=1.0$')
#plt.legend()
if save:
    plt.savefig('figures/Ediff_g1_p2.pdf')

# ------- Iteration Convergence ------- #
plt.figure(figsize=[7,5])
plt.title('Convergence with alpha (Iterations)',fontsize=20)
plt.xlabel('Number of basis states',fontsize=16)
plt.ylabel('Number of iterations',fontsize=16)
#plt.plot(x,c02,label=r'CCD $\alpha=0.2$')
plt.plot(x,c04,label=r'CCD $\alpha=0.4$')
#plt.plot(x,c06,label=r'CCD $\alpha=0.6$')
plt.plot(x,c08,label=r'CCD $\alpha=0.8$')
plt.plot(x,c10,label=r'CCD $\alpha=1.0$')
plt.legend()
if save:
    plt.savefig('figures/conv_g1_p2.pdf')


# ---------- Time ----------- #
plt.figure(figsize=[7,5])
plt.title('Time',fontsize=20)
plt.xlabel('Number of basis states',fontsize=16)
plt.ylabel('CPU time',fontsize=16)
#plt.semilogy(x,t02,label=r'CCD $\alpha=0.2$')
plt.semilogy(x,t04,label=r'CCD $\alpha=0.4$')
#plt.semilogy(x,t06,label=r'CCD $\alpha=0.6$')
plt.semilogy(x,t08,label=r'CCD $\alpha=0.8$')
plt.semilogy(x,t10,label=r'CCD $\alpha=1.0$')
plt.semilogy(x,CI_t_g1,label='CI')
plt.legend()
if save:
    plt.savefig('figures/time_g1_p2.pdf')

plt.show()
