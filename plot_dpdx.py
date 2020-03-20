
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#############################################
# Read in data from an ASCII data table
# WRITE (1,'(I6,E17.8,9E22.13)') j, epp, dedx_a_e, dedx_a_i, dedx_a_tot, &
# dedxc_a_e, dedxc_a_i, dedxc_a_tot, dedxq_a_e, dedxq_a_i, dedxq_a_tot
#
data1 = np.genfromtxt('dpdx_1.out')
E = data1[:, 1]          
dedxe = data1[:, 2]
dedxi = data1[:, 3]
dedx = data1[:, 4]
dedxce = data1[:, 5]
dedxci = data1[:, 6]
dedxc = data1[:, 7]
dedxqe = data1[:, 8]
dedxqi = data1[:, 9]
dedxq = data1[:, 10]

# plot stopping power total
plt.plot(E,dedxe, label=r'$dE_e/dx$')
plt.plot(E,dedxi, label=r'$dE_I/dx$')
plt.plot(E,dedx, label=r'$dE_{tot}/dx$')
#
xmax = 3.5e3
ymax = 5.0e7
plt.xlim(0, xmax)
plt.ylim(0, ymax)
plt.xlabel(r'$E \,\, {\rm [KeV]}$')
plt.ylabel(r'$dE/dx \,\, {\rm [KeV/cm]}$')
plt.title(r'$dE_b/dx$: BPS')
plt.legend(loc=5)
plt.grid(True)
plt.savefig('plot_dpdx_1.png')
plt.show()

#############################################
# Read in data from an ASCII data table
# WRITE (2,'(I6,E17.8,9E22.13)') j, epp, c_e, c_i, c_tot, cc_e, cc_i, cc_tot, &
# cq_e, cq_i, cq_tot           
#
data2 = np.genfromtxt('dpdx_2.out')
E = data2[:, 1]          
Clle = data2[:, 2]
Clli = data2[:, 3]
Cll = data2[:, 4]
Cllce = data2[:, 5]
Cllci = data2[:, 6]
Cllc = data2[:, 7]
Cllqe = data2[:, 8]
Cllqi = data2[:, 9]
Cllq = data2[:, 10]

# plot stopping power total
plt.plot(E,Clle, label='dedx_e')
plt.plot(E,Clli, label='dedx_i')
plt.plot(E,Cll, label='dedx_tot')
#
xmax = 3.5e3
ymax = 5.0e7
plt.xlim(0, xmax)
#plt.ylim(0, ymax)
plt.xlabel(r'$E \,\, {\rm [KeV]}$')
plt.ylabel(r'$C^{\ell\ell} \,\, {\rm [keV/c/cm]}$')
plt.title(r'$C^{\ell\ell}_b$: BPS')
plt.legend(loc=5)
plt.grid(True)
plt.savefig('plot_dpdx_2.png')
plt.show()


#############################################
# Read in data from an ASCII data table
# WRITE (2,'(I6,E17.8,9E22.13)') j, epp, c_e, c_i, c_tot, cc_e, cc_i, cc_tot, &
# cq_e, cq_i, cq_tot           
#
data3 = np.genfromtxt('dpdx_3.out')
E2 = data3[:, 1]          
vpdPdxe2 = data3[:, 2]
vpdPdxi2 = data3[:, 3]
vpdPdx2 = data3[:, 4]
vpdPdxce2 = data3[:, 5]
vpdPdxci2 = data3[:, 6]
vpdPdxc2 = data3[:, 7]
vpdPdxqe2 = data3[:, 8]
vpdPdxqi2 = data3[:, 9]
vpdPdxq2 = data3[:, 10]

# plot stopping power total
plt.plot(E2,vpdPdxe2, label=r'${\bf v} \cdot d{\bf P}_e/dx$')
plt.plot(E2,vpdPdxi2, label=r'${\bf v} \cdot d{\bf P}_I/dx$')
plt.plot(E2,vpdPdx2, label=r'${\bf v} \cdot d{\bf P}_{tot}/dx$')
#
xmax = 3.5e3
ymax = 5.0e7
plt.xlim(0, xmax)
#plt.ylim(0, ymax)
plt.xlabel(r'$E \,\, {\rm [KeV]}$')
plt.ylabel(r'${\bf v}p \cdot d{\bf P}_b/dx \,\, {\rm [keV^2/s/cm^2]}$')
plt.title(r'${\bf v} \cdot d{\bf P}_b/dx$: BPS')
plt.legend(loc=5)
plt.grid(True)
plt.savefig('plot_dpdx_3.png')
plt.show()


