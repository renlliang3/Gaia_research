import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
import astropy.coordinates as asc
import numpy as np
import random
import csv
import os

script_dir = os.path.dirname(__file__)
rel_path = "Data/A_corrected_mean_AG_EBminR-vrad.fits"
abs_file_path = os.path.join(script_dir, rel_path)

readresults = Table.read(abs_file_path,format='fits')
results = np.array(readresults)
distances = 1000/results['parallax']

#for i in range(len(distances)):
#	if np.isnan(results['radial_velocity'][i]) == True:
#		results['radial_velocity'][i] = 0

coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='icrs', obstime='J2015.5')

coordinates_galactic = coordinates_ICRS.galactic
#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, frame='galactic', obstime='J2015.5')

coordinates_cartesian = np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))

x = coordinates_cartesian[:,0]
y = coordinates_cartesian[:,1]
z = coordinates_cartesian[:,2]

#distance = np.sqrt(x**2+y**2+z**2)

#indices650 = [i for i in range(len(distance)) if distance[i] <= 650]

#x = np.array([x[i] for i in indices650])
#y = np.array([y[i] for i in indices650])
#z = np.array([z[i] for i in indices650])

vx = coordinates_galactic.velocity.d_x.value
vy = coordinates_galactic.velocity.d_y.value
vz = coordinates_galactic.velocity.d_z.value

towrite = np.column_stack((x,y,z,vx,vy,vz)).tolist()
towrite.insert(0,['xg','yg','zg', 'vx','vy','vz'])

with open('A_corrected_mean_AG_EBminR-vrad_v_xyz.csv', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
    wr.writerows(towrite)
