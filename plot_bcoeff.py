import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#############################################
# Read in data from an ASCII data table
data = np.genfromtxt('bcoeff_1.out')
E = data[:, 1]     
Be = data[:, 2]    
Bi = data[:, 3]    
B = data[:, 4]     
Be_c = data[:, 5]  
Bi_c = data[:, 6]  
B_c = data[:, 7]   
Be_q = data[:, 8]  
Bi_q = data[:, 9]  
B_q = data[:, 10]  

# plot stopping power total
plt.plot(E,B, label='B_tot')
plt.plot(E,Be, label='B_e')
plt.plot(E,Bi, label='B_i')

#
xmax = 3.5e3
ymax = 0.1e10
ymin =-0.6e10
plt.xlim(0, xmax)
plt.ylim(ymin, ymax)
plt.xlabel(r'$E \,\, {\rm [keV]}$')
plt.ylabel(r'$B \,\, {\rm [keV \cdot s/cm^2]}$')
plt.title(r'$B$: BPS')
plt.legend(loc=4)
plt.grid(True)
plt.savefig('plot_bcoeff_1.png')
plt.show()


#############################################
# plot stopping power classical
plt.plot(E,B_c, label='B_tot_c')
plt.plot(E,Be_c, label='B_e_c')
plt.plot(E,Bi_c, label='B_i_c')

#
xmax = 3.5e3
ymax = 0.1e10
ymin =-0.6e10
plt.xlim(0, xmax)
plt.ylim(ymin, ymax)
plt.xlabel(r'$E \,\, {\rm [keV]}$')
plt.ylabel(r'$B^C \,\, {\rm [keV \cdot s/cm^2]}$')
plt.title(r'classical: $B^C}$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_bcoeff_2.png')
plt.show()


#############################################
# plot stopping power quantum
plt.plot(E,B_q, label='B_tot_q')
plt.plot(E,Be_q, label='B_e_q')
plt.plot(E,Bi_q, label='B_i_q')

#
xmax = 3.5e3
ymax = 1.0e6
ymin = 0
plt.xlim(0, xmax)
plt.ylim(ymin, ymax)
plt.xlabel(r'$E \,\, {\rm [keV]}$')
plt.ylabel(r'$B^Q \,\, {\rm [keV \cdot s/cm^2]}$')
plt.title(r'quantum: $B^{Q}$: BPS')
plt.legend(loc=5)
plt.grid(True)
plt.savefig('plot_bcoeff_3.png')
plt.show()


#############################################
# Read in data from an ASCII data table
data = np.genfromtxt('bcoeff_2.out') # dE_perp/dx
E = data[:, 1]     
dEperpdxe = data[:, 2]    
dEperpdxi = data[:, 3]    
dEperpdx = data[:, 4]     
dEperpdxe_c = data[:, 5]  
dEperpdxi_c = data[:, 6]  
dEperpdx_c = data[:, 7]   
dEperpdxe_q = data[:, 8]  
dEperpdxi_q = data[:, 9]  
dEperpdx_q = data[:, 10]  

# plot stopping power total
plt.plot(E,dEperpdx, label='dEperpdx_tot')
plt.plot(E,dEperpdxe, label='dEperpdx_e')
plt.plot(E,dEperpdxi, label='dEperpdx_i')

#
xmax = 3.5e3
#ymax = 0.1e10
#ymin =-0.6e10
plt.xlim(0, xmax)
#plt.ylim(ymin, ymax)
plt.xlabel(r'$E \,\, {\rm [keV]}$')
plt.ylabel(r'$dEperpdx \,\, {\rm [keV/cm]}$')
plt.title(r'$dEperpdx$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_bcoeff_4.png')
plt.show()
