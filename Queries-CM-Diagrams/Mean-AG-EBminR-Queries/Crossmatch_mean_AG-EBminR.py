import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.io.votable import parse
from astropy.table import Table
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

extra_data = np.genfromtxt("Data/Mean_extinction_excess_SourceID_500pc_AG-4.4_BminR-1.7.dat", names = True, dtype=None)
source_id = extra_data['source_id']

readresults = Table.read('Data/All-star-500pc-AG-4.4-BminR-1.7.fits',format='fits')
results = np.array(readresults)

nonmatches = np.array([])
j=0

for i in range(len(source_id)):
	not_found = True
	while not_found:
		if source_id[i]!=results['source_id'][j]:
			nonmatches = np.append(nonmatches,j)
			j+=1
		else:
			not_found = False
			j+=1
	print(j)

for k in range(len(nonmatches)):
	results = np.delete(results,nonmatches[k]-k)

j=0
nonmatches_check = np.array([])
for i in range(len(source_id)):
	not_found = True
	while not_found:
		if source_id[i]!=results['source_id'][j]:
			nonmatches_check = np.append(nonmatches_check,j)
			j+=1
		else:
			not_found = False
			j+=1
	print(j)
