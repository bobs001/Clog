
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#############################################
# Read in data from an ASCII data table
data = np.genfromtxt('ccoeff_1.out')
E = data[:, 1]       
Clle = data[:, 2]    
Clli = data[:, 3]    
Cll = data[:, 4]     
Clle_c = data[:, 5]  
Clli_c = data[:, 6]  
Cll_c = data[:, 7]   
Clle_q = data[:, 8]  
Clli_q = data[:, 9]  
Cll_q = data[:, 10]  

#############################################
# plot stopping power total
plt.plot(E,Clle, label='C_e')
plt.plot(E,Clli, label='C_i')
plt.plot(E,Cll, label='C_tot')

#
xmax = 3.5e3
ymax = 5.0
plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [keV]}$')
plt.ylabel(r'$C^{\ell\ell} \,\, {\rm [keV^2 \cdot s/cm^2]}$')
plt.title(r'$C^{\ell\ell}$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_1.png')
plt.show()


#############################################
# plot stopping power classical
plt.plot(E,Clle_c, label='C_e_c')
plt.plot(E,Clli_c, label='C_i_c')
plt.plot(E,Cll_c, label='C_tot_c')

#
xmax = 3.5e3
ymax = 5.0
plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [keV]}$')
plt.ylabel(r'$C^{\ell\ell} \,\, {\rm [keV^2 \cdot s/cm^2]}$')
plt.title(r'classical: $C^{\ell\ell, C}$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_2.png')
plt.show()

#############################################
# plot stopping power quantum
plt.plot(E,Clle_q, label='C_e_q')
plt.plot(E,Clli_q, label='C_i_q')
#plt.plot(E,Cll_q, label='C_tot_q')

#
xmax = 3.5e3
ymax = 5.0
plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [keV]}$')
plt.ylabel(r'$C^{\ell\ell} \,\, {\rm [keV^2 \cdot s/cm^2]}$')
plt.title(r'quantum: $C^{\ell\ell, Q}$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_3.png')
plt.show()
