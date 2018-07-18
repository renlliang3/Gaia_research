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

"""
readresults = Table.read('Data/A-OB-vrad-500pc-mean_AG_EBminR_v_xyz.fits',format='fits')
results = np.array(readresults)

x, y, z, vx, vy, vz = results['xg'], results['yg'], results['zg'], results['vx']+14.0, results['vy']+12.24, results['vz']+7.25

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
"""

readresults_OD = Table.read('Velacontour4.5e10-7.dat',format='ascii')
results_OD = np.array(readresults_OD)

x_OD, y_OD, z_OD, vx_OD, vy_OD, vz_OD = results_OD['xg'], results_OD['yg'], results_OD['zg'], results_OD['vx']+14.0, results_OD['vy']+12.24, results_OD['vz']+7.25

v_xyz_OD = np.vstack([vx_OD,vy_OD,vz_OD])

bw_OD = 3

kde_fit_OD = KernelDensity(bandwidth=bw_OD, kernel='gaussian').fit(v_xyz_OD.T)

rangesize_OD=100
boxes_OD=100j

VX_OD, VY_OD, VZ_OD = np.mgrid[-rangesize_OD:rangesize_OD:100j, -rangesize_OD:rangesize_OD:100j, -rangesize_OD:rangesize_OD:100j]
positions_OD = np.vstack([VX_OD.ravel(), VY_OD.ravel(), VZ_OD.ravel()])
dens_OD = np.reshape(np.exp(kde_fit_OD.score_samples(positions_OD.T)), VX_OD.shape)

#Visualize the density estimate as isosurfaces
figure = mlab.figure('myfig')
figure.scene.disable_render = True # Super duper trick
#mlab.contour3d(VX, VY, VZ, dens, contours=[0.25*10**-5,0.5*10**-5,1*10**-5,2*10**-5,2.5*10**-5,3*10**-5], opacity=0.3, colormap = 'Blues', figure = figure)
mlab.contour3d(VX_OD, VY_OD, VZ_OD, dens_OD, contours=[0.1*10**-4,0.5*10**-4,1*10**-4,1.5*10**-4,2*10**-4,2.5*10**-4], opacity=0.3, vmin=np.min(dens_OD), vmax=np.max(dens_OD), colormap = 'Reds', figure = figure)
mlab.axes(xlabel='U (km/s)', ylabel='V (km/s)', zlabel='W (km/s)', nb_labels = 3, figure = figure)
#mlab.quiver3d(x, y, z, vx, vy, vz)

#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)
mlab.points3d(vx_OD, vy_OD, vz_OD, np.full(len(vx_OD),1), scale_factor = 1, color=(0,0,0), figure = figure)

figure.scene.disable_render = False # Super duper trick
mlab.show()
