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

readresults = Table.read('Data/Updated-A-Query-Radial-Velocities_v_xyz.fits',format='fits')
results = np.array(readresults)

readresults2 = Table.read('Sourceswithincontour.dat',format='ascii')
results2 = np.array(readresults2)

x = results['xg']
y = results['yg']
z = results['zg']

vx = results['vx']
vy = results['vy']
vz = results['vz']

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
mlab.contour3d(VX, VY, VZ, dens, contours=[2*10**-5,4*10**-5,6*10**-5], opacity=0.3, colormap = 'Blues', figure = figure)
#mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes(xlabel='vx (km/s)', ylabel='vy (km/s)', zlabel='vz (km/s)', nb_labels = 5, figure = figure)
#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)

figure.scene.disable_render = False # Super duper trick
mlab.show()
