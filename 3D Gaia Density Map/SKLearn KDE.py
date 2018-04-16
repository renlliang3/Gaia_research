import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from sklearn.neighbors import KernelDensity
from mayavi import mlab
import astropy.coordinates as asc
import matplotlib.pyplot as plt
import numpy as np

readresults = Table.read('Gaiadr1.gaia_source with min parallax 1.54 and error over parallax less than 0.2 and M_G less than 2.vot',format='votable')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, frame='icrs', obstime='J2015')
coordinates_cartesian = np.column_stack((coordinates_ICRS.cartesian.x.value, coordinates_ICRS.cartesian.y.value, coordinates_ICRS.cartesian.z.value))
#coordinates_galactic = coordinates_ICRS.galactic
#coordinates_cartesian_galactic_converted = np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))

#Check if their galactic coordinates matches ours 
#realgalactic = np.column_stack((results['l'],results['b'],distances))
#convertedgalactic = np.column_stack((coordinates_galactic.l.degree,coordinates_galactic.b.degree,coordinates_galactic.distance.pc))
#differencegalactic = convertedgalactic-realgalactic

#Check for differences in the cartesian coordinates between using the converted and the real galactic coordinates
#coordinates_cartesian_galactic_real = asc.spherical_to_cartesian(distances*u.pc, results['b']*u.degree, results['l']*u.degree)
#coordinates_cartesian_galactic_real = np.column_stack((coordinates_cartesian_galactic_real[0].value, coordinates_cartesian_galactic_real[1].value, coordinates_cartesian_galactic_real[2].value))
#differencecartesian = coordinates_cartesian_galactic_real - coordinates_cartesian_galactic_converted

x = coordinates_cartesian[:,0]#.filled(0)
y = coordinates_cartesian[:,1]#.filled(0)
z = coordinates_cartesian[:,2]#.filled(0)

xyz = np.vstack([x,y,z])

d = xyz.shape[0]
n = xyz.shape[1]
bw = (n*(d+2)/4.)**(-1./(d+4))

#kde_fit = KernelDensity(bandwidth=bw, kernel='epanechnikov').fit(xyz.T)
kde_fit = KernelDensity(3, kernel='gaussian').fit(xyz.T)

X, Y, Z=np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]
positions=np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
#dens doesn't converge look up how to make the shapes of the numpy arrays fit by
#checking the example script that did work
#dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), X.shape)

# Visualize the density estimate as isosurfaces
#mlab.contour3d(X, Y, Z, dens, opacity=0.5)
#mlab.axes()
#mlab.show()

#plt.plot
#plt.imshow(np.rot90(dens), cmap=plt.cm.viridis, extent=[x.min(), x.max(), y.min(), y.max(), z.min(), z.max()], aspect='auto')
#plt.contour(X_gaia, Y_gaia, dens_gaia, colors='k', linewidth=0.01, levels = [10,30,50])
#plt.text(-0.1, 10.0, 'Gaia DR1')
#plt.xlim(-0.3,2.0)
#plt.ylim(-4.0,12.5)
#plt.xlabel('(B-V)')
#plt.yticks([])
#plt.gca().invert_yaxis()

#plt.savefig('KDE')

#fig = plt.figure()
#ax = Axes3D(fig)

#ax.scatter(coordinates_cartesian_galactic_real[:,0], coordinates_cartesian_galactic_real[:,1], coordinates_cartesian_galactic_real[:,2], s=0.1)
#plt.show()
