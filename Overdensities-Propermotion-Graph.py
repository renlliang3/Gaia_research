import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from matplotlib.colors import LogNorm
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
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, frame='icrs', obstime='J2015.5')
#coordinates_ICRS = asc.ICRS(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr)
coordinates_galactic = coordinates_ICRS.galactic

#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='galactic', obstime='J2015.5')

#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, frame='galactic', obstime='J2015.5')

coordinates_cartesian = np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))

x = coordinates_cartesian[:,0]#.filled(0)
y = coordinates_cartesian[:,1]#.filled(0)
z = coordinates_cartesian[:,2]#.filled(0)

#counts,xbins,ybins,image = plt.hist2d(distances*(4.74 * 10**-3)*results['pmra'],distances*(4.74 * 10**-3)*results['pmdec'],bins=60,normed=True,norm=LogNorm(), cmap = 'Blues')
#counts,xbins,ybins,image = plt.hist2d(results['pmra'],results['pmdec'],bins=40,normed=True,norm=LogNorm(),cmap = 'Blues')
hb = plt.hexbin(distances*(4.74 * 10**-3)*results['pmra'], distances*(4.74 * 10**-3)*results['pmdec'], extent=(-50,50,-50,50), gridsize=80, bins='log', cmap = 'Blues')
#hb = plt.hexbin(results['pmra'], results['pmdec'], gridsize=80, extent=(-50,50,-50,50), bins='log',  cmap = 'Blues')
plt.colorbar()
#plt.contour(counts.transpose(), extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()], colors='k', linewidth=0.01), levels = [0.001])
#plt.text(-0.1, 10.0, 'Gaia DR1')
plt.xlim(-50,50)
plt.ylim(-50,50)
plt.xlabel(r'$V_{Tra} \ (km/s)$')
plt.ylabel(r'$V_{Tdec} \ (km/s)$')
#plt.xlabel(r'$\mu_{ra} \ (mas/yr)$')
#plt.ylabel(r'$\mu_{dec} \ (mas/yr)$')
#plt.savefig('Proper-Motion-Katz-contour-0.1.png')
plt.savefig('Tangential-Velocites-Katz-Contour-0.1.png')
