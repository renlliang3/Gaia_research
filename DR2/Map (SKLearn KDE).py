import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from sklearn.neighbors import KernelDensity
from mayavi import mlab
import astropy.coordinates as asc
import numpy as np
import random

readresults = Table.read('Reduced_DR2_Radial_Velocity.csv',format='csv')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='icrs', obstime='J2015.5')
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

vx = coordinates_ICRS.velocity.d_x.value
vy = coordinates_ICRS.velocity.d_y.value
vz = coordinates_ICRS.velocity.d_z.value

random_indices = random.sample(range(len(x)),5000)

x = np.array([x[i] for i in sorted(random_indices)])
y = np.array([y[i] for i in sorted(random_indices)])
z = np.array([z[i] for i in sorted(random_indices)])

vx = np.array([vx[i] for i in sorted(random_indices)])
vy = np.array([vy[i] for i in sorted(random_indices)])
vz = np.array([vz[i] for i in sorted(random_indices)])

xyz = np.vstack([x,y,z])

d = xyz.shape[0]
n = xyz.shape[1]
#bw = (n*(d+2)/4.)**(-1./(d+4))
bw = 28

kde_fit = KernelDensity(bandwidth=bw, kernel='gaussian').fit(xyz.T)

X, Y, Z = np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]
positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), X.shape)

#Visualize the density estimate as isosurfaces
mlab.contour3d(X, Y, Z, dens, opacity=0.5)
mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes()
mlab.show()
