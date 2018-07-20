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

readresults = Table.read("Data/A-OB-500pc-maxvtan40-mean_AG_EBminR.fits",format='fits')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, frame='icrs', obstime='J2015.5')
#coordinates_ICRS = asc.ICRS(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr)
coordinates_galactic = coordinates_ICRS.galactic

pm_b = coordinates_galactic.pm_b.value
pm_l_cosb = coordinates_galactic.pm_l_cosb.value

vtb = distances*(4.74057 * 10**-3)*pm_b
vtl = distances*(4.74057 * 10**-3)*pm_l_cosb

indices_low_vtan = np.where((vtb<=10)&(vtb>=-20))

readresults = readresults[indices_low_vtan]

readresults.write('A-OB-500pc-maxvtan40-10vb-20-mean_AG_EBminR.fits', format='fits')
