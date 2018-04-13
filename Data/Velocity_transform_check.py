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

readresults = Table.read('GOG_GUMS_DR2_7e6.fits',format='fits')
results = np.array(readresults)

dist = results['dist']
list = []

for i in range(len(dist)):
	if (dist[i] < 650):
		list.append(i)

list = np.array(list)

dist = np.array([dist[i] for i in list])

source_id = np.array([results['source_id'][i] for i in list])
ra = np.array([results['ra'][i] for i in list])
dec = np.array([results['dec'][i] for i in list])
varpi = np.array([results['varpi'][i] for i in list])
e_varpi = np.array([results['e_varpi'][i] for i in list])
pmra = np.array([results['pmra'][i] for i in list])
e_pmra = np.array([results['e_pmra'][i] for i in list])
pmdec = np.array([results['pmdec'][i] for i in list])
e_pmdec = np.array([results['e_pmdec'][i] for i in list])
vrad = ((60/22)**0.5)*np.array([results['vrad'][i] for i in list])
e_vrad = np.array([results['e_vrad'][i] for i in list])
G_mag = np.array([results['G_Mag'][i] for i in list])
GRVS_Mag = np.array([results['GRVS_Mag'][i] for i in list])
GBP_Mag = np.array([results['GBP_Mag'][i] for i in list])
GRP_Mag = np.array([results['GRP_Mag'][i] for i in list])
dist = np.array([results['dist'][i] for i in list])

pmra_cosdec = pmra*np.cos(dec*np.pi/180)
                           
coordinates_ICRS = asc.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, distance=dist*u.pc, pm_ra_cosdec=pmra_cosdec*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vrad*u.km/u.s, frame='icrs', obstime='J2010')

x = coordinates_ICRS.cartesian.x.value#.filled(0)
y = coordinates_ICRS.cartesian.y.value#.filled(0)
z = coordinates_ICRS.cartesian.z.value#.filled(0)

vx = coordinates_ICRS.velocity.d_x.value
vy = coordinates_ICRS.velocity.d_y.value
vz = coordinates_ICRS.velocity.d_z.value

v_Tra = pmra_cosdec * dist * 4.74 * 10**-3 #mas/yr -> km/s
v_Tdec = pmdec * dist * 4.74 * 10**-3 #mas/yr -> km/s

vx_check = (vrad * np.cos(dec*np.pi/180) * np.cos(ra*np.pi/180)) - (v_Tra * np.sin(ra*np.pi/180)) - (v_Tdec * np.sin(dec*np.pi/180) * np.cos(ra*np.pi/180))
vy_check = (vrad * np.cos(dec*np.pi/180) * np.sin(ra*np.pi/180)) + (v_Tra * np.cos(ra*np.pi/180)) - (v_Tdec * np.sin(dec*np.pi/180) * np.sin(ra*np.pi/180))
vz_check = (vrad * np.sin(dec*np.pi/180)) + (v_Tdec * np.cos(dec*np.pi/180))

random_indices = random.sample(range(len(x)),100000)

x = np.array([x[i] for i in sorted(random_indices)])
y = np.array([y[i] for i in sorted(random_indices)])
z = np.array([z[i] for i in sorted(random_indices)])

vx = np.array([vx[i] for i in sorted(random_indices)])
vy = np.array([vy[i] for i in sorted(random_indices)])
vz = np.array([vz[i] for i in sorted(random_indices)])

print('Mean of (vx^2+vy^2+vz^2)^0.5: '+str(np.mean((vx**2+vy**2+vz**2)**0.5)))
print('Mean of (vx_check^2+vy_check^2+vz_check^2)^0.5: '+str(np.mean((vx_check**2+vy_check**2+vz_check**2)**0.5)))
print('SD of (vx^2+vy^2+vz^2)^0.5: '+str(np.var((vx**2+vy**2+vz**2)**0.5)**0.5))
print('SD of (vx_check^2+vy_check^2+vz_check^2)^0.5: '+str(np.var((vx_check**2+vy_check**2+vz_check**2)**0.5)**0.5))

print('Mean of vx: '+str(np.mean(vx)))
print('Mean of vx_check: '+str(np.mean(vx_check)))
print('sd of vx: '+str(np.var(vx)**0.5))
print('sd of vx_check: '+str(np.var(vx_check)**0.5))
