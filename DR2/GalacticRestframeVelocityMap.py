import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.coordinates import LSR
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from mayavi import mlab
import astropy.coordinates as asc
import matplotlib.pyplot as plt
import numpy as np
import random
import csv

readresults = Table.read('Reduced_DR2_Radial_Velocity.csv',format='csv')
results = np.array(readresults)

distances = 1000/results['parallax']
list650 = []
list1000 = []

for i in range(len(distances)):
	if (distances[i] < 650) & (results['parallax_over_error'][i] >= 5.0):
		list650.append(i)
	if (distances[i] < 1000) & (results['parallax_over_error'][i] >= 5.0):
		list1000.append(i)

list = np.array(list650)

source_id = np.array([results['source_id'][i] for i in list])
ra = np.array([results['ra'][i] for i in list])
ra_error = np.array([results['ra_error'][i] for i in list])
dec = np.array([results['dec'][i] for i in list])
dec_error = np.array([results['dec_error'][i] for i in list])
parallax = np.array([results['parallax'][i] for i in list])
parallax_error = np.array([results['parallax_error'][i] for i in list])
parallax_over_error = np.array([results['parallax_over_error'][i] for i in list])

pmra = np.array([results['pmra'][i] for i in list])
pmra_error = np.array([results['pmra_error'][i] for i in list])
pmdec = np.array([results['pmdec'][i] for i in list])
pmdec_error = np.array([results['pmdec_error'][i] for i in list])

radial_velocity = np.array([results['radial_velocity'][i] for i in list])
radial_velocity_error = np.array([results['radial_velocity_error'][i] for i in list])

g_mean_mag = np.array([results['phot_g_mean_mag'][i] for i in list])
bp_mean_mag = np.array([results['phot_bp_mean_mag'][i] for i in list])
rp_mean_mag = np.array([results['phot_rp_mean_mag'][i] for i in list])

distances = np.array([distances[i] for i in list])
                           
coordinates_ICRS = asc.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, distance=distances*u.pc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=radial_velocity*u.km/u.s, frame='icrs', obstime='J2015.5')

gc1 = coordinates_ICRS.transform_to(asc.Galactocentric)

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

mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes()
mlab.show()

#x = coordinates_ICRS.cartesian.x.value#.filled(0)
#y = coordinates_ICRS.cartesian.y.value#.filled(0)
#z = coordinates_ICRS.cartesian.z.value#.filled(0)

#vx = coordinates_ICRS.velocity.d_x.value
#vy = coordinates_ICRS.velocity.d_y.value
#vz = coordinates_ICRS.velocity.d_z.value

#v_Tra = pmra_cosdec * dist * 4.74 * 10**-3 #mas/yr -> km/s
#v_Tdec = pmdec * dist * 4.74 * 10**-3 #mas/yr -> km/s

#vx_check = (vrad * np.cos(dec) * np.cos(ra)) - (v_Tra * np.sin(ra)) - (v_Tdec * np.sin(dec) * np.cos(ra))
#vy_check = (vrad * np.cos(dec) * np.sin(ra)) + (v_Tra * np.cos(ra)) - (v_Tdec * np.sin(dec) * np.sin(ra))
#vz_check = (vrad * np.sin(dec)) + (v_Tdec * np.cos(dec))

#random_indices = random.sample(range(len(x)),100000)

#x = np.array([x[i] for i in sorted(random_indices)])
#y = np.array([y[i] for i in sorted(random_indices)])
#z = np.array([z[i] for i in sorted(random_indices)])

#vx = np.array([vx[i] for i in sorted(random_indices)])
#vy = np.array([vy[i] for i in sorted(random_indices)])
#vz = np.array([vz[i] for i in sorted(random_indices)])

#mlab.quiver3d(x, y, z, vx, vy, vz)
#mlab.axes()
#mlab.show()

#xyz = np.vstack([x,y,z])

#kde = stats.gaussian_kde(xyz)

#X, Y, Z = np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]
