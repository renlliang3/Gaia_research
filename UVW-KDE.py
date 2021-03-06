import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from sklearn.neighbors import KernelDensity
from mayavi import mlab
import astropy.coordinates as asc
import numpy as np
import random
import os

readresults = Table.read('Data/A-OB-vrad-500pc-mean_AG_EBminR_v_xyz.fits',format='fits')
results = np.array(readresults)

readresults2 = Table.read('Orioncontour1e10-8.dat',format='ascii')
results2 = np.array(readresults2)

x, y, z, vx, vy, vz = results['xg'], results['yg'], results['zg'], results['vx']+14.0, results['vy']+12.24, results['vz']+7.25

x_OD, y_OD, z_OD, vx_OD, vy_OD, vz_OD = results2['xg'], results2['yg'], results2['zg'], results2['vx']+14.0, results2['vy']+12.24, results2['vz']+7.25

#random_indices = random.sample(range(len(x)),5000)

#x = np.array([x[i] for i in sorted(random_indices)])
#y = np.array([y[i] for i in sorted(random_indices)])
#z = np.array([z[i] for i in sorted(random_indices)])

#vx = np.array([vx[i] for i in sorted(random_indices)])
#vy = np.array([vy[i] for i in sorted(random_indices)])
#vz = np.array([vz[i] for i in sorted(random_indices)])

v_xyz = np.vstack([vx,vy,vz])

d = v_xyz.shape[0]
n = v_xyz.shape[1]
#bw = (n*(d+2)/4.)**(-1./(d+4))
bw = 3

kde_fit = KernelDensity(bandwidth=bw, kernel='gaussian').fit(v_xyz.T)

rangesize=100
boxes=100j
VX, VY, VZ = np.mgrid[-rangesize:rangesize:boxes, -rangesize:rangesize:boxes, -rangesize:rangesize:boxes]
positions = np.vstack([VX.ravel(), VY.ravel(), VZ.ravel()])
dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), VX.shape)

binsize = 2*rangesize/(boxes.imag)
w = int(rangesize/binsize)

#Visualize the density estimate as isosurfaces
figure = mlab.figure('myfig')
figure.scene.disable_render = True # Super duper trick
mlab.contour3d(VX, VY, VZ, dens, contours=[0.25*10**-5,0.5*10**-5,1*10**-5,2*10**-5,2.5*10**-5,3*10**-5], opacity=0.3, colormap = 'Blues', figure = figure)
#mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes(xlabel='U (km/s)', ylabel='V (km/s)', zlabel='W (km/s)', nb_labels = 5, figure = figure)
#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)
mlab.points3d(vx_OD, vy_OD, vz_OD, np.full(len(vx_OD),1), scale_factor = 1, color=(0,0,0), figure = figure)

figure.scene.disable_render = False # Super duper trick
mlab.show()
