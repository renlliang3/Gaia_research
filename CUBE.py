import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from astropy.io import fits
from  numpy.lib.recfunctions import append_fields
#from astropy.io import ascii
#import imageio.core
import astropy.coordinates as asc
import numpy as np
import random
import csv
import os
import sys

#Source_ID
readresults = Table.read("Data/All-star-500pc-AG-4.4-BminR-1.7.fits",format='fits')
results = np.array(readresults)

#mean_AG
fp_hdu_AG = fits.open("Data/mean_AG_map.fits")
fp_AG = fp_hdu_AG[0].data

#mean_AG
fp_hdu_EBminR = fits.open("Data/mean_EBminR_map.fits")
fp_EBminR = fp_hdu_EBminR[0].data

binSize = 10
w = 500//binSize
z_height = 500//binSize

#galactic coordinates
fp_hdu3 = fits.open("Data/All-star-500pc-AG-4.4-BminR-1.7_xyz.fits")
fp3 = fp_hdu3[1].data

xgal, ygal, zgal = fp3['xg'], fp3['yg'], fp3['zg']

#Code
source_id_arr = np.array([],dtype='i8')
AG_arr = np.array([],dtype='f8')
EBminR_arr = np.array([],dtype='f8')
#Check in which box every star falls then append to a new array the source_id of that star and the corresponding mean_AG and mean_EBminR
print('Making Map')
N = len(xgal)
for  i in range(N):
	x = int(round(float(xgal[i])/binSize))+w
	y = int(round(float(ygal[i])/binSize))+w
	z = int(round(float(zgal[i])/binSize))+z_height
	if x >= 0 and x < 2*w and y >= 0 and y < 2*w and z > 0 and z < 2*z_height:
		source_id_arr = np.append(source_id_arr,results['source_id'][i])
		AG_arr = np.append(AG_arr,fp_AG[x][y][z])
		EBminR_arr = np.append(EBminR_arr,fp_EBminR[x][y][z])
		print(i)	

source_id_arr = np.rec.array(source_id_arr, dtype=[('source_id', np.int64)])
new_arr = append_fields(source_id_arr, 'AG', AG_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'EBminR', EBminR_arr, usemask=False, dtypes=[np.float64])

np.savetxt('Mean_extinction_excess_SourceID_500pc_AG-4.4_BminR-1.7.dat', new_arr, header='source_id AG EBminR', fmt = '%2li %1.8f %1.8f')
