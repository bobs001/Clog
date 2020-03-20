
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

#############################################
# Read in data from an ASCII data table
data = np.genfromtxt('bcoeff_1.out')
E = data[:, 1]          # [keV]
Be = data[:, 2]      # [MeV/mu-m]
Bi = data[:, 3]      # [MeV/mu-m]
B = data[:, 4]       # [MeV/mu-m]
Be_c = data[:, 5]      # [MeV/mu-m]
Bi_c = data[:, 6]      # [MeV/mu-m]
B_c = data[:, 7]       # [MeV/mu-m]
Be_q = data[:, 8]      # [MeV/mu-m]
Bi_q = data[:, 9]      # [MeV/mu-m]
B_q = data[:, 10]       # [MeV/mu-m]

#############################################
# plot stopping power total
plt.plot(E,Be, label='B_e')
plt.plot(E,Bi, label='B_i')
plt.plot(E,B, label='B_tot')

xmax = 3.5
ymax = 5.0
#plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$B \,\, {\rm [MeV/\mu m]}$*')
plt.title(r'$B$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_1.png')
plt.show()


#############################################
# plot stopping power classical
plt.plot(E,Be_c, label='B_e_c')
plt.plot(E,Bi_c, label='B_i_c')
plt.plot(E,B_c, label='B_tot_c')

#
xmax = 3.5
ymax = 5.0
#plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$B^{C} \,\, {\rm [MeV/\mu m]}$*')
plt.title(r'classical: $B^C}$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_2.png')
plt.show()

#############################################
# plot stopping power quantum
plt.plot(E,Be_q, label='B_e_q')
plt.plot(E,Bi_q, label='B_i_q')
#plt.plot(E,B_q, label='B_tot_q')

#
xmax = 3.5
ymax = 5.0
#plt.xlim(0, xmax)
#plt.ylim(-8, 5)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$B^{Q} \,\, {\rm [MeV^2/\mu m]}$**')
plt.title(r'quantum: $B^{Q}$: BPS')
plt.legend(loc=0)
plt.grid(True)
plt.savefig('plot_ccoeff_3.png')
plt.show()
