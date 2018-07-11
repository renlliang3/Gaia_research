import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from mayavi import mlab
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter as gf
from scipy.special import expit
import matplotlib.pyplot as plt
import astropy.coordinates as asc
import numpy as np
import random
import csv
import os

contour = np.genfromtxt("Data/OB-Katz-contour-0.1.dat", names = True, dtype=None)
contour_id = contour['source_id']

readresults = Table.read("Data/OB-Katz.fits",format='fits')
results = np.array(readresults)

matches = np.array([])
j=0
for i in range(len(contour_id)):
	not_found = True
	while not_found:
		if contour_id[i]==results['source_id'][j]:
			matches = np.append(matches,j)
			not_found = False
		else:
			j+=1

newresults = np.append(results[int(matches[0])],results[int(matches[1])])

k = 2
while k <= len(matches)-1:
	newresults = np.append(newresults,results[int(matches[k])])
	k+=1

results = newresults
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
#coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, frame='icrs', obstime='J2015.5')
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, frame='icrs', obstime='J2015.5')
coordinates_galactic = coordinates_ICRS.galactic

#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='galactic', obstime='J2015.5')

#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, frame='galactic', obstime='J2015.5')

coordinates_cartesian = np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))

xgal = coordinates_cartesian[:,0]#.filled(0)
ygal = coordinates_cartesian[:,1]#.filled(0)
zgal = coordinates_cartesian[:,2]#.filled(0)

binSize = 3
w = 500//binSize
z_height = 350//binSize

### Create map with the bright and young stars

sigma = 5
spread = 50

a = np.zeros((w*2,w*2,2*z_height))
print('Making Map')
N = len(xgal)
for  i in range(N):
    x = int(round(float(xgal[i])/binSize))+w
    y = int(round(float(ygal[i])/binSize))+w
    z = int(round(float(zgal[i])/binSize))+z_height
    if x >= 0 and x < 2*w and y >= 0 and y < 2*w and z >= 0 and z < 2*z_height:
        a[x][y][z] += 1

print('Smoothing')
gaussian = gf(a, sigma=sigma, truncate=3)
b = 2*(expit(spread*gaussian)-0.5)

from mayavi import mlab
import  moviepy.editor as mpy

fig_myv = mlab.figure(size=(750,700), bgcolor=(1, 1, 1), fgcolor = (0,0,0))
#duration = 9
mlab.clf()
src = mlab.pipeline.scalar_field(b)
m = mlab.pipeline.iso_surface(src,
                              #contours = [0.2, 0.3, 0.5],
                              contours = [0.1],
                              opacity = 0.3, colormap = 'Reds', figure = fig_myv)
mlab.axes(#ranges = [-1300, 1300, -1300, 1300, -3500, 3500],
          ranges = [-500, 500, -500, 500, -350, 350],
          color = (0,0,0), nb_labels = 5, xlabel = 'X(pc)', ylabel = 'Y(pc)', zlabel = 'Z(pc)', figure = fig_myv)
fig_myv.scene.disable_render = False
mlab.view(azimuth=0, elevation= 0,figure=fig_myv)
mlab.show()
