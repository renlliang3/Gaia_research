import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.coordinates import GalacticLSR
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

readresults = Table.read("Data/A-OB-500pc-maxvtan40-10vb-20-mean_AG_EBminR.fits",format='fits')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, frame='icrs', obstime='J2015.5')
#coordinates_ICRS = asc.ICRS(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr)
coordinates_galactic = coordinates_ICRS.galactic

pm_dec = results['pmdec']
pm_ra_cosdec = results['pmra']

vtra = distances*(4.74057 * 10**-3)*pm_ra_cosdec
vtdec = distances*(4.74057 * 10**-3)*pm_dec

vt_radec_total = np.sqrt(vtra**2+vtdec**2)

pm_b = coordinates_galactic.pm_b.value
pm_l_cosb = coordinates_galactic.pm_l_cosb.value

vtb = distances*(4.74057 * 10**-3)*pm_b
vtl = distances*(4.74057 * 10**-3)*pm_l_cosb

def densKDE2D(x, y, bw, limits, boxes):
	xy = np.vstack([x,y])
	kde_fit_xy = KernelDensity(bandwidth=bw, kernel='gaussian').fit(xy.T)
	
	X, Y = np.mgrid[limits[0]:limits[1]:boxes, limits[2]:limits[3]:boxes]
	positions = np.vstack([X.ravel(), Y.ravel()])
	return np.reshape(np.exp(kde_fit_xy.score_samples(positions.T)), X.shape)

dens = densKDE2D(vtl, vtb, 0.25, [-40,40,-20,10], 200j)

#counts,xbins,ybins,image = plt.hist2d(vtl,vtb,bins=100,Norm = colors.PowerNorm(0.5), cmap = 'Blues')
#counts,xbins,ybins,image = plt.hist2d(results['pmra'],results['pmdec'],bins=40,normed=True,norm=LogNorm(),cmap = 'Blues')
#hb = plt.hexbin(vtb, vtl, gridsize=600, Norm = colors.PowerNorm(0.5), bins='log', cmap = 'Blues')
#hb = plt.hexbin(vtb, vtl, gridsize=600, Norm = colors.PowerNorm(0.5), bins='log', cmap = 'Blues')
#hb = plt.hexbin(results['pmra'], results['pmdec'], gridsize=80, extent=(-50,50,-50,50), bins='log',  cmap = 'Blues')

plt.imshow(np.rot90(dens), norm = colors.PowerNorm(0.5), cmap='Blues', extent = [-40, 40, -20, 10])
#plt.contour(counts.transpose(), extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()], colors='k', linewidth=0.01), levels = [0.001])
#plt.text(-0.1, 10.0, 'Gaia DR1')
plt.colorbar()
#plt.xlim(-100,100)
#plt.ylim(-100,100)
plt.xlabel(r'$V_{T_l} \ (km/s)$')
plt.ylabel(r'$V_{T_b} \ (km/s)$')
#plt.xlabel(r'$\mu_{ra} \ (mas/yr)$')
#plt.ylabel(r'$\mu_{dec} \ (mas/yr)$')
#plt.savefig('Proper-Motion-Katz-contour-0.1.png')
plt.savefig('Tangential-Galactic-Velocites-KDE-imshow-200boxes-A-OB-500pc-maxvtan40-10vb-20-mean_AG_EBminR.png')
