import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from sklearn.neighbors import KernelDensity
from matplotlib.colors import LogNorm
from astropy.io import ascii
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
ascii.write(results, 'OB-Katz-contour-0.1.dat')
