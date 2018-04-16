import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from mayavi import mlab
import astropy.coordinates as asc
import matplotlib.pyplot as plt
import numpy as np
import random

readresults = Table.read('Gaiadr1.gaia_source with min parallax 1.54 and error over parallax less than 0.2 and M_G less than 2.vot',format='votable')
results = np.array(readresults)
distances = 1000/results['parallax']
pmra_cosdec = results['pmra']*np.cos(results['dec']*np.pi/180)

dec = results['dec']
dec_error = results['dec_error']
ra = results['ra']
ra_error = results['ra_error']

indices = []
for i in range(len(dec)):
    if (np.abs(dec_error[i]/dec[i]) < 0.1 and np.abs(ra_error[i]/ra[i]) < 0.1):
            indices.append(i)

newdec = np.array([dec[i] for i in indices])
newdec_error = np.array([dec_error[i] for i in indices])
newra = np.array([ra[i] for i in indices])
newra_error = np.array([ra_error[i] for i in indices])

plt.subplot(2,2,1)
plt.plot(dec_error/dec)
plt.subplot(2,2,2)
plt.plot(ra_error/ra)
plt.subplot(2,2,3)
plt.plot(newdec_error/newdec)
plt.subplot(2,2,4)
plt.plot(newra_error/newra)
plt.show()
