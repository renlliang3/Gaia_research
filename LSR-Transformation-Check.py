import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.coordinates import LSR
from astropy.coordinates import GalacticLSR
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
import astropy.coordinates as asc
import numpy as np
import random
import csv
import os

readresults = Table.read("Data/A-OB-vrad-500pc-mean_AG_EBminR.fits",format='fits')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian.
#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, frame='galactic')
#coordinates_cartesian_check= np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))

#x_check = coordinates_cartesian_check[:,0]
#y_check = coordinates_cartesian_check[:,1]
#z_check = coordinates_cartesian_check[:,2]

#converted_ra = coordinates_galactic.icrs.ra.value
#converted_dec = coordinates_galactic.icrs.dec.value

#coordinates_ICRS_and_vel = asc.SkyCoord(ra=converted_ra*u.degree, dec=converted_dec*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='icrs', obstime='J2015.5')

coordinates_ICRS_and_vel = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='icrs', obstime='J2015.5')

coordinates_galactic_and_vel = coordinates_ICRS_and_vel.galactic
#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, frame='galactic', obstime='J2015.5')

coordinates_cartesian = np.column_stack((coordinates_galactic_and_vel.cartesian.x.value, coordinates_galactic_and_vel.cartesian.y.value, coordinates_galactic_and_vel.cartesian.z.value))

source_id = results['source_id']
x = coordinates_cartesian[:,0]
y = coordinates_cartesian[:,1]
z = coordinates_cartesian[:,2]

#distance = np.sqrt(x**2+y**2+z**2)

#indices650 = [i for i in range(len(distance)) if distance[i] <= 650]

#x = np.array([x[i] for i in indices650])
#y = np.array([y[i] for i in indices650])
#z = np.array([z[i] for i in indices650])

vx = coordinates_galactic_and_vel.velocity.d_x.value
vy = coordinates_galactic_and_vel.velocity.d_y.value
vz = coordinates_galactic_and_vel.velocity.d_z.value

vx_LSR_check = vx + 11.1
vy_LSR_check = vy + 12.24
vz_LSR_check = vz + 7.25

coordinates_GalacticLSR = coordinates_ICRS_and_vel.transform_to(GalacticLSR)

vx_LSR = coordinates_GalacticLSR.velocity.d_x.value
vy_LSR = coordinates_GalacticLSR.velocity.d_y.value
vz_LSR = coordinates_GalacticLSR.velocity.d_z.value
