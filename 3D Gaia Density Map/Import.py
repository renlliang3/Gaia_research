import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from mpl_toolkits.mplot3d import Axes3D
import astropy.coordinates as asc
import matplotlib.pyplot as plt
import numpy as np

readresults = Table.read('Gaiadr1.gaia_source with min parallax 1.54 and error over parallax less than 0.2 and M_G less than 1.7.vot',format='votable')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, frame='icrs', obstime='J2015')
coordinates_galactic = coordinates_ICRS.galactic
coordinates_cartesian_galactic_converted = np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))

#Check if their galactic coordinates matches ours 
realgalactic = np.column_stack((results['l'],results['b'],distances))
convertedgalactic = np.column_stack((coordinates_galactic.l.degree,coordinates_galactic.b.degree,coordinates_galactic.distance.pc))
differencegalactic = convertedgalactic-realgalactic

#Check for differences in the cartesian coordinates between using the converted and the real galactic coordinates
coordinates_cartesian_galactic_real = asc.spherical_to_cartesian(distances*u.pc, results['b']*u.degree, results['l']*u.degree)
coordinates_cartesian_galactic_real = np.column_stack((coordinates_cartesian_galactic_real[0].value, coordinates_cartesian_galactic_real[1].value, coordinates_cartesian_galactic_real[2].value))
differencecartesian = coordinates_cartesian_galactic_real - coordinates_cartesian_galactic_converted

fig = plt.figure()
ax = Axes3D(fig)

ax.scatter(coordinates_cartesian_galactic_real[0:20,0], coordinates_cartesian_galactic_real[0:20,1], coordinates_cartesian_galactic_real[0:20,2])
plt.show()
