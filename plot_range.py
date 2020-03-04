import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams

# WRITE (6,'(I6,E17.8,6D22.13)') j, epp, dedx_tot, dedx_i, dedx_e
# Read in data from an ASCII data table
data = np.genfromtxt('dedx.001.out')
E = data[:, 1]/1000.    # [MeV]
dedxe = data[:, 2]      # [MeV/mu-m]
dedxi = data[:, 3]      # [MeV/mu-m]
dedx = data[:, 4]       # [MeV/mu-m
# Create a loglog plot of data
plt.plot(E,dedxe)
plt.plot(E,dedxi)
plt.plot(E,dedx)
xmax = 3.5
ymax = 0.3
plt.xlim(0, xmax)
plt.ylim(0, ymax)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Stopping Power: BPS')
plt.legend(loc=2)
plt.grid(True)
plt.savefig('plot_range_1.png')
plt.show()

#                                0  1      2      3      4   5         6       7
# WRITE (6,'(I6,E17.8,6D22.13)') j, tt(j), xt(j), vt(j), ep, dedx_tot, dedx_i, dedx_e
data = np.genfromtxt('range.001.out')
x = data[:, 2]*10.E4 # mu-m
dedx = data[:, 5]    # MeV/mu-m
dedxi = data[:, 6]   # MeV/mu-m
dedxe = data[:, 7]   # MeV/mu-m
# # Create a loglog plot of data
plt.plot(x, dedxe)
plt.plot(x, dedxi)
plt.plot(x, dedx)
xmax = 0.0032
xmax = 320.0 # mu-m
plt.xlim(0, xmax)
plt.ylim(0, ymax)
plt.xlabel(r'$x \,\, {\rm [cm]}$')
plt.ylabel(r'$dE/dx \,\, {\rm [MeV/\mu m]}$')
plt.title('Stopping Power: BPS')
plt.legend(loc=2)
plt.grid(True)
plt.savefig('plot_range_2.png')
plt.show()
