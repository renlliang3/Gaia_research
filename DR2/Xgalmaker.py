import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from mayavi import mlab
import astropy.coordinates as asc
import numpy as np
import random
import csv

readresults = Table.read('Reduced_DR2_Radial_Velocity.csv',format='csv')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='icrs', obstime='J2015.5')

coordinates_cartesian = np.column_stack((coordinates_ICRS.cartesian.x.value, coordinates_ICRS.cartesian.y.value, coordinates_ICRS.cartesian.z.value))

x = coordinates_cartesian[:,0]#.filled(0)
y = coordinates_cartesian[:,1]#.filled(0)
z = coordinates_cartesian[:,2]#.filled(0)

vx = coordinates_ICRS.velocity.d_x.value
vy = coordinates_ICRS.velocity.d_y.value
vz = coordinates_ICRS.velocity.d_z.value


towrite = np.column_stack((x,y,z,vx,vy,vz)).tolist()
towrite.insert(0,['xg','yg','zg','vx','vy','vz'])

with open('Reduced_DR2_Radial_Velocity_xyz', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
    wr.writerows(towrite)
