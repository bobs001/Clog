
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#############################################
# Read in data from an ASCII data table
data1 = np.genfromtxt('acoeff.dedx_1.out')
E1 = data1[:, 1]          # [keV]
dedxe1 = data1[:, 2]      # [MeV/mu-m]
dedxi1 = data1[:, 3]      # [MeV/mu-m]
dedx1 = data1[:, 4]       # [MeV/mu-m]
dedxce1 = data1[:, 5]      # [MeV/mu-m]
dedxci1 = data1[:, 6]      # [MeV/mu-m]
dedxc1 = data1[:, 7]       # [MeV/mu-m]
dedxqe1 = data1[:, 8]      # [MeV/mu-m]
dedxqi1 = data1[:, 9]      # [MeV/mu-m]
dedxq1 = data1[:, 10]       # [MeV/mu-m]
#
data2= np.genfromtxt('acoeff.dedx_2.out')
E2 = data2[:, 1]          # [keV]
dedxe2 = data2[:, 2]      # [MeV/mu-m]
dedxi2 = data2[:, 3]      # [MeV/mu-m]
dedx2 = data2[:, 4]       # [MeV/mu-m]
dedxce2 = data2[:, 5]      # [MeV/mu-m]
dedxci2 = data2[:, 6]      # [MeV/mu-m]
dedxc2 = data2[:, 7]       # [MeV/mu-m]
dedxqe2 = data2[:, 8]      # [MeV/mu-m]
dedxqi2 = data2[:, 9]      # [MeV/mu-m]
dedxq2 = data2[:, 10]       # [MeV/mu-m]

data3= np.genfromtxt('acoeff.dedx_3.out')
E3 = data3[:, 1]          # [keV]
del_dedxe3 = data3[:, 2]      # [MeV/mu-m]
del_dedxi3 = data3[:, 3]      # [MeV/mu-m]
del_dedx3 = data3[:, 4]       # [MeV/mu-m]

#############################################
# plot stopping power total
plt.plot(E1,dedxe1, label='dedx_e')
plt.plot(E1,dedxi1, label='dedx_i')
plt.plot(E1,dedx1, label='dedx_tot')
plt.plot(E2,dedxe2)
plt.plot(E2,dedxi2)
plt.plot(E2,dedx2)
#
xmax = 3.5
ymax = 5.0
plt.xlim(0, xmax)
plt.ylim(0, ymax)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Stopping Power: BPS')
plt.legend(loc=2)
plt.grid(True)
plt.savefig('plot_acoeff.dedx_1.png')
plt.show()

#############################################
# plot stopping power classical
#
plt.plot(E1,dedxce1, label='dedxc_e')
plt.plot(E1,dedxci1, label='dedxc_i')
plt.plot(E1,dedxc1, label='dedxc_tot')
plt.plot(E2,dedxce2)
plt.plot(E2,dedxci2)
plt.plot(E2,dedxc2)
#
xmax = 3.5
ymax = 8.0
plt.xlim(0, xmax)
plt.ylim(0, ymax)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Stopping Power: BPS Classical')
plt.legend(loc=2)
plt.grid(True)
plt.savefig('plot_acoeff.dedx_2.png')
plt.show()

#############################################
# plot stopping power quantum
#
plt.plot(E1,dedxqe1, label='dedxq_e')
plt.plot(E1,dedxqi1, label='dedxq_i')
plt.plot(E1,dedxq1, label='dedxq_tot')
plt.plot(E2,dedxqe2)
plt.plot(E2,dedxqi2)
plt.plot(E2,dedxq2)
#
xmax = 3.5
ymin = -4.0
ymax = 0
plt.xlim(0, xmax)
plt.ylim(ymin, ymax)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Stopping Power: BPS Quantum')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_acoeff.dedx_3.png')
plt.show()

#############################################
# 
plt.plot(E1, del_dedxe3, label='del_dedx_e')
plt.plot(E1, del_dedxi3, label='del_dedx_i')
plt.plot(E1, del_dedx3, label='del_dedx_tot')
xmax = 3.5
#ymin =-0.1e-3
plt.xlim(0, xmax)
#plt.ylim(-4.e-3, 4.e-3)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'error')
plt.title('Error: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_acoeff.dedx_4.png')
plt.show()

##### is thie the right A-coeff?
#############################################
# Read in data from an ASCII data table
data = np.genfromtxt('acoeff_dedx_0.out')
E = data[:, 1]          # [keV]
Ae = data[:, 2]      # [MeV/mu-m]
Ai = data[:, 3]      # [MeV/mu-m]
A_tot = data[:, 4]       # [MeV/mu-m]
Ae_c = data[:, 5]      # [MeV/mu-m]
Ai_c = data[:, 6]      # [MeV/mu-m]
A_c = data[:, 7]       # [MeV/mu-m]
Ae_q = data[:, 8]      # [MeV/mu-m]
Ai_q = data[:, 9]      # [MeV/mu-m]
A_q = data[:, 10]       # [MeV/mu-m]

#############################################
# plot stopping power total
plt.plot(E,Ae, label='A_e')
plt.plot(E,Ai, label='A_i')
plt.plot(E,A_tot, label='A_tot')

xmax = 3.5
ymax = 5.0
plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$C^{\ell \ell} \,\, {\rm [MeV/\mu m]}$')
plt.title(r'${\calA}$ profile: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_5.png')
plt.show()
