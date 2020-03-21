
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#############################################
# Read in data from an ASCII data table
data = np.genfromtxt('pcoeff_1.out')
E = data[:, 1]       
Plle = data[:, 2]    
Plli = data[:, 3]    
Pll = data[:, 4]     
Plle_c = data[:, 5]  
Plli_c = data[:, 6]  
Pll_c = data[:, 7]   
Plle_q = data[:, 8]  
Plli_q = data[:, 9]  
Pll_q = data[:, 10]  

# plot stopping power total
plt.plot(E,Plle, label='P_e')
plt.plot(E,Plli, label='P_i')
plt.plot(E,Pll, label='P_tot')

#
xmax = 3.5e3
ymax = 5.0e7
plt.xlim(0, xmax)
plt.ylim(0, ymax)
plt.xlabel(r'$E \,\, {\rm [keV]}$')
plt.ylabel(r'$v^k \, dP^k/dx \,\, {\rm [keV/cm]}$')
plt.title(r'$v^k dP^k/dx$: BPS')
plt.legend(loc=5)
plt.grid(True)
plt.savefig('plot_pcoeff_1.png')
plt.show()

