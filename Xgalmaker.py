import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from astropy.io import fits
import astropy.coordinates as asc
import numpy as np
import random
import csv
import os

script_dir = os.path.dirname(__file__)
rel_path = "Data/A-OB-500pc-maxvtan40-10vb-20-mean_AG_EBminR.fits"
abs_file_path = os.path.join(script_dir, rel_path)

readresults = Table.read(abs_file_path,format='fits')
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

coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, frame='icrs', obstime='J2015.5')

#coordinates_ICRS_and_vel = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='icrs', obstime='J2015.5')

#coordinates_galactic_and_vel = coordinates_ICRS_and_vel.galactic
coordinates_galactic = coordinates_ICRS.galactic
#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, frame='galactic', obstime='J2015.5')

#coordinates_cartesian = np.column_stack((coordinates_galactic_and_vel.cartesian.x.value, coordinates_galactic_and_vel.cartesian.y.value, coordinates_galactic_and_vel.cartesian.z.value))
coordinates_cartesian = np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))

source_id = results['source_id']
x = coordinates_cartesian[:,0]
y = coordinates_cartesian[:,1]
z = coordinates_cartesian[:,2]

pm_b = coordinates_galactic.pm_b.value
pm_l_cosb = coordinates_galactic.pm_l_cosb.value

vtb = distances*(4.74057 * 10**-3)*pm_b
vtl = distances*(4.74057 * 10**-3)*pm_l_cosb

#distance = np.sqrt(x**2+y**2+z**2)

#indices650 = [i for i in range(len(distance)) if distance[i] <= 650]

#x = np.array([x[i] for i in indices650])
#y = np.array([y[i] for i in indices650])
#z = np.array([z[i] for i in indices650])

#vx = coordinates_galactic_and_vel.velocity.d_x.value
#vy = coordinates_galactic_and_vel.velocity.d_y.value
#vz = coordinates_galactic_and_vel.velocity.d_z.value

#towrite = np.column_stack((source_id,x,y,z)).tolist()#,vx,vy,vz)).tolist()
#towrite = list(zip(source_id.tolist(),x.tolist(),y.tolist(),z.tolist(),vx.tolist(),vy.tolist(),vz.tolist()))
#towrite.insert(0,['source_id','xg','yg','zg','vx','vy','vz'])

towrite = list(zip(source_id.tolist(),x.tolist(),y.tolist(),z.tolist(),distances.tolist(),pm_l_cosb.tolist(),pm_b.tolist(),vtl.tolist(),vtb.tolist()))
towrite.insert(0,['source_id','xg','yg','zg','distances','pm_l_cosb','pm_b','vtl','vtb'])

with open('A-OB-500pc-maxvtan40-10vb-20-mean_AG_EBminR_xyz_pm_vtan.csv', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
    wr.writerows(towrite)
