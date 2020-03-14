
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#############################################
# Read in data from an ASCII data table
data = np.genfromtxt('pcoeff_1.out')
E = data[:, 1]          # [keV]
Plle = data[:, 2]      # [MeV/mu-m]
Plli = data[:, 3]      # [MeV/mu-m]
Pll = data[:, 4]       # [MeV/mu-m]
Plle_c = data[:, 5]      # [MeV/mu-m]
Plli_c = data[:, 6]      # [MeV/mu-m]
Pll_c = data[:, 7]       # [MeV/mu-m]
Plle_q = data[:, 8]      # [MeV/mu-m]
Plli_q = data[:, 9]      # [MeV/mu-m]
Pll_q = data[:, 10]       # [MeV/mu-m]

#############################################
# plot stopping power total
plt.plot(E,Plle, label='P_e')
plt.plot(E,Plli, label='P_i')
plt.plot(E,Pll, label='P_tot')

#
xmax = 3.5
ymax = 5.0
plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$P^{\ell \ell} \,\, {\rm [MeV/\mu m]}$*')
plt.title(r'$P^{\ell\ell}$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_1.png')
plt.show()

