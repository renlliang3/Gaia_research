import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from mayavi import mlab
import astropy.coordinates as asc
import numpy as np
import random

readresults = Table.read('GOG_GUMS_DR2_7e6.fits',format='fits')
results = np.array(readresults)

dist = results['dist']
list = []

for i in range(len(dist)):
	if dist[i]<650:
		list.append(i)

list = np.array(list)

dist = [dist[i] for i in list]

source_id = [results['source_id'][i] for i in list]
ra = [results['ra'][i] for i in list]
dec = np.array([results['dec'][i] for i in list])
varpi = [results['varpi'][i] for i in list]
e_varpi = [results['e_varpi'][i] for i in list]
pmra = np.array([results['pmra'][i] for i in list])
e_pmra = [results['e_pmra'][i] for i in list]
pmdec = [results['pmdec'][i] for i in list]
e_pmdec = [results['e_pmdec'][i] for i in list]
vrad = [results['vrad'][i] for i in list]
e_vrad = [results['e_vrad'][i] for i in list]
G_mag = [results['G_Mag'][i] for i in list]
GRVS_Mag = [results['GRVS_Mag'][i] for i in list]
GBP_Mag = [results['GBP_Mag'][i] for i in list]
GRP_Mag = [results['GRP_Mag'][i] for i in list]
dist = [results['dist'][i] for i in list]

pmra_cosdec = pmra*np.cos(dec*np.pi/180)

coordinates_ICRS = asc.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, distance=dist*u.pc, pm_ra_cosdec=pmra_cosdec*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vrad*u.km/u.s, frame='icrs', obstime='J2010')
#coordinates_cartesian = np.column_stack((coordinates_ICRS.cartesian.x.value, coordinates_ICRS.cartesian.y.value, coordinates_ICRS.cartesian.z.value))
gc1 = coordinates_ICRS.transform_to(asc.Galactocentric)

#x = coordinates_cartesian[:,0]#.filled(0)
#y = coordinates_cartesian[:,1]#.filled(0)
#z = coordinates_cartesian[:,2]#.filled(0)

#vx = coordinates_ICRS.velocity.d_x.value
#vy = coordinates_ICRS.velocity.d_y.value
#vz = coordinates_ICRS.velocity.d_z.value

x = gc1.cartesian.x.value#.filled(0)
y = gc1.cartesian.y.value#.filled(0)
z = gc1.cartesian.z.value#.filled(0)

vx = gc1.v_x.value
vy = gc1.v_y.value
vz = gc1.v_z.value

random_indices = random.sample(range(len(x)),1000)

x = np.array([x[i] for i in sorted(random_indices)])
y = np.array([y[i] for i in sorted(random_indices)])
z = np.array([z[i] for i in sorted(random_indices)])

vx = np.array([vx[i] for i in sorted(random_indices)])
vy = np.array([vy[i] for i in sorted(random_indices)])
vz = np.array([vz[i] for i in sorted(random_indices)])

print("1")

xyz = np.vstack([x,y,z])

kde = stats.gaussian_kde(xyz)

print("2")

X, Y, Z = np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]

coords = np.vstack([item.ravel() for item in [X, Y, Z]])
density = kde(coords).reshape(X.shape)

print("3")

# Visualize the density estimate as isosurfaces
mlab.contour3d(X, Y, Z, density, opacity=0.5)
mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes()
mlab.show()
