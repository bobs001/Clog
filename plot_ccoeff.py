
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#############################################
# Read in data from an ASCII data table
data = np.genfromtxt('ccoeff_dedx_1.out')
E = data[:, 1]          # [keV]
Clle = data[:, 2]      # [MeV/mu-m]
Clli = data[:, 3]      # [MeV/mu-m]
Cll = data[:, 4]       # [MeV/mu-m]
Clle_c = data[:, 5]      # [MeV/mu-m]
Clli_c = data[:, 6]      # [MeV/mu-m]
Cll_c = data[:, 7]       # [MeV/mu-m]
Clle_q = data[:, 8]      # [MeV/mu-m]
Clli_q = data[:, 9]      # [MeV/mu-m]
Cll_q = data[:, 10]       # [MeV/mu-m]

#############################################
# plot stopping power total
plt.plot(E,Clle, label='C_e')
plt.plot(E,Clli, label='C_i')
plt.plot(E,Cll, label='C_tot')

#
xmax = 3.5
ymax = 5.0
plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$C^{\ell \ell} \,\, {\rm [MeV/\mu m]}$*')
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
xmax = 3.5
ymax = 5.0
plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$C^{\ell \ell} \,\, {\rm [MeV/\mu m]}$*')
plt.title(r'$C^{\ell\ell}$: BPS')
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
xmax = 3.5
ymax = 5.0
plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$C^{\ell \ell} \,\, {\rm [MeV/\mu m]}$*')
plt.title(r'$C^{\ell\ell}$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_3.png')
plt.show()
