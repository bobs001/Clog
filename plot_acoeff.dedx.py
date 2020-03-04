
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc, rcParams


# Read in data from an ASCII data table
data1 = np.genfromtxt('acoeff.dedx_1.out')
E1 = data1[:, 1]          # [keV]
dedxe1 = data1[:, 2]      # [MeV/mu-m]
dedxi1 = data1[:, 3]      # [MeV/mu-m]
dedx1 = data1[:, 4]       # [MeV/mu-m]
#
data2= np.genfromtxt('acoeff.dedx_2.out')
E2 = data2[:, 1]          # [keV]
dedxe2 = data2[:, 2]      # [MeV/mu-m]
dedxi2 = data2[:, 3]      # [MeV/mu-m]
dedx2 = data2[:, 4]       # [MeV/mu-m]


# Create a loglog plot of data
plt.plot(E1,dedxe1)
plt.plot(E1,dedxi1)
plt.plot(E1,dedx1)
plt.plot(E2,dedxe2)
plt.plot(E2,dedxi2)
plt.plot(E2,dedx2)

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
plt.clf()

# error
#
E3 = (E1 - E2) / E1
dedxe3 = (dedxe1 - dedxe2 ) / dedxe1
dedxi3 = (dedxi1 - dedxi2 ) / dedxi1
dedx3 = (dedx1 - dedx2 ) / dedx1

# Create a loglog plot of data
plt.plot(E3,dedxe3)
plt.plot(E3,dedxi3)
plt.plot(E3,dedx3)

print "***", dedx3

xmax = 3.5
ymin =-0.1e-3
plt.xlim(0, xmax)
#plt.ylim(-4.e-3, 4.e-3)
plt.xlabel(r'$E \,\, {\rm [MeV]}$')
plt.ylabel(r'error')
plt.title('Stopping Power: BPS')
plt.legend(loc=2)
plt.grid(True)
#plt.savefig('plot_acoeff.dedx_2.png')
plt.show()
