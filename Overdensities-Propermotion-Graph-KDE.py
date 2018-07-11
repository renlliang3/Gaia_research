import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from sklearn.neighbors import KernelDensity
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
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

vtra = distances*(4.74 * 10**-3)*results['pmra']
vtdec = distances*(4.74 * 10**-3)*results['pmdec']

vtradec = np.vstack([vtra,vtdec])
#pmradec = np.vstack([results['pmra'],results['pmdec']])

d = vtradec.shape[0]
n = vtradec.shape[1]
#bw = (n*(d+2)/4.)**(-1./(d+4))
bw = 0.25
#kde_fit = KernelDensity(bandwidth=bw, kernel='gaussian').fit(pmradec.T)
kde_fit_vtan = KernelDensity(bandwidth=bw, kernel='gaussian').fit(vtradec.T)

#PMRA, PMDEC = np.mgrid[results['pmra'].min():results['pmra'].max():1000j, results['pmdec'].min():results['pmdec'].max():1000j]
#positions = np.vstack([PMRA.ravel(), PMDEC.ravel()])
#dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), PMRA.shape)

VTRA, VTDEC = np.mgrid[vtra.min():vtra.max():500j, vtdec.min():vtdec.max():500j]
positions = np.vstack([VTRA.ravel(), VTDEC.ravel()])
dens = np.reshape(np.exp(kde_fit_vtan.score_samples(positions.T)), VTRA.shape)

#plt.imshow(np.rot90(dens), norm=colors.LogNorm(vmin=10**-7, vmax=dens.max()), cmap='Blues', extent=[results['pmra'].min(),results['pmra'].max(),results['pmdec'].min(),results['pmdec'].max()])
plt.imshow(np.rot90(dens), norm=colors.LogNorm(vmin=10**-6, vmax=dens.max()), cmap='Blues', extent=[vtra.min(),vtra.max(),vtdec.min(),vtdec.max()])
plt.colorbar()

plt.xlim(-50,50)
plt.ylim(-50,50)
plt.xlabel(r'$V_{Tra} \ (km/s)$')
plt.ylabel(r'$V_{Tdec} \ (km/s)$')

#plt.xlabel(r'$\mu_{ra} \ (mas/yr)$')
#plt.ylabel(r'$\mu_{dec} \ (mas/yr)$')
#plt.savefig('Proper-Motion-Katz-contour-0.3-KDE-test.png')
plt.savefig('Tangential-Velocites-Katz-Contour-0.1-KDE-bw=0.25-boxes=500-range50.png')
