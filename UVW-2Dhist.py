import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from sklearn.neighbors import KernelDensity
from mayavi import mlab
from matplotlib.colors import LogNorm
import astropy.coordinates as asc
import matplotlib.pyplot as plt
import numpy as np
import random
import os

script_dir = os.path.dirname(__file__)
rel_path = "Data/A_corrected_mean_AG_EBminR-vrad_v_xyz.fits"
abs_file_path = os.path.join(script_dir, rel_path)

readresults = Table.read(abs_file_path,format='fits')
results = np.array(readresults)

x = results['xg']
y = results['yg']
z = results['zg']

vx = results['vx']
vy = results['vy']
vz = results['vz']

counts_gaia,xbins_gaia,ybins_gaia,image_gaia = plt.hist2d(vx,vy,bins=60,norm=LogNorm(), cmap = 'Blues')
plt.colorbar()
#plt.contour(counts_gaia.transpose(), extent=[xbins_gaia.min(),xbins_gaia.max(),ybins_gaia.min(),ybins_gaia.max()], colors='k', linewidth=0.01, levels = [0.1])
#plt.text(-0.1, 10.0, 'Gaia DR1')
plt.xlim(-100,100)
plt.ylim(-100,100)
plt.xlabel(r'$v_x$')
plt.ylabel(r'$v_y$')

plt.savefig('2dhist')
