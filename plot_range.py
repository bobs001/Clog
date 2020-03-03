import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams


# Read in data from an ASCII data table
data = np.genfromtxt('dedx.001.out')
E = data[:, 1]/1000.    # [MeV]
dedxe = data[:, 2]      # [MeV/mu-m]
dedxi = data[:, 3]      # [MeV/mu-m]
dedx = data[:, 4]       # [MeV/mu-m
print "*** dE/dx **", dedxe
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
