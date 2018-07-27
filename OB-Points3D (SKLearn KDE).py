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

# OB ass positions
ob = np.genfromtxt("Data/OB_positions.dat", names = True, dtype=None)
x_OB, y_OB, z_OB = ob['xg'], ob['yg'], ob['zg']
OB_Names = ob['OB_Names']

readresults = Table.read("Data/A-OB-500pc-maxvtan40-10vb-20-mean_AG_EBminR_xyz_pm_vtan_vrad0.fits",format='fits')
results = np.array(readresults)

#x, y, z, vx, vy, vz = results['xg'], results['yg'], results['zg'], results['vx']+14.0, results['vy']+12.24, results['vz']+7.25
x, y, z, vx, vy, vz = results['xg'], results['yg'], results['zg'], results['vx'], results['vy'], results['vz']

xyz = np.vstack([x[0:5000],y[0:5000],z[0:5000]])

d = xyz.shape[0]
n = xyz.shape[1]
#bw = (n*(d+2)/4.)**(-1./(d+4))
bw = 20

kde_fit = KernelDensity(bandwidth=bw, kernel='gaussian').fit(xyz.T)

crange=500

X, Y, Z = np.mgrid[-crange:crange:100j, -crange:crange:100j, -crange:crange:100j]
positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), X.shape)

"""
readresults_OD = Table.read('Sourceswithincontour3e10-5.dat',format='ascii')
results_OD = np.array(readresults_OD)

x_OD, y_OD, z_OD, vx_OD, vy_OD, vz_OD = results_OD['xg'], results_OD['yg'], results_OD['zg'], results_OD['vx']+14.0, results_OD['vy']+12.24, results_OD['vz']+7.25

xyz_OD = np.vstack([x_OD,y_OD,z_OD])

bw_OD = 20

kde_fit_OD = KernelDensity(bandwidth=bw_OD, kernel='gaussian').fit(xyz_OD.T)

crange_OD = 500

X_OD, Y_OD, Z_OD = np.mgrid[-crange_OD:crange_OD:100j, -crange_OD:crange_OD:100j, -crange_OD:crange_OD:100j]
positions_OD = np.vstack([X_OD.ravel(), Y_OD.ravel(), Z_OD.ravel()])
dens_OD = np.reshape(np.exp(kde_fit_OD.score_samples(positions_OD.T)), X_OD.shape)
"""

#Visualize the density estimate as isosurfaces
figure = mlab.figure('myfig')
#figure.scene.disable_render = True # Super duper trick
#mlab.contour3d(X, Y, Z, dens, extent = [-650, 650, -650, 650, -650, 650], opacity=0.3, colormap = 'Blues', figure = figure)
mlab.contour3d(X, Y, Z, dens, contours=[7.45*10**-9,7.5*10*--9,10*10**-9,12.5*10**-9,15*10**-9,17.5*10**-9,20*10**-9,22.5*10**-9], opacity=0.3, vmin=np.min(dens), vmax=np.max(dens), colormap = 'Blues', figure = figure)
#mlab.contour3d(X_OD, Y_OD, Z_OD, dens_OD, contours=[0.5*10**-7,1.5*10**-7,2.5*10**-7,3.5*10**-7,4.5*10**-7], opacity=0.3, vmin=np.min(dens_OD), vmax=np.max(dens_OD), colormap = 'Reds', figure = figure)
mlab.quiver3d(x[0:1000], y[0:1000], z[0:1000], vx[0:1000]+14.0, vy[0:1000]+12.24, vz[0:1000]+7.25)
#mlab.quiver3d(x[0:5000], y[0:5000], z[0:5000], vx[0:5000], vy[0:5000], vz[0:5000])
mlab.axes(nb_labels = 3, figure = figure)
#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)
mlab.points3d(x_OB, y_OB, z_OB, np.full(len(x_OB),10), scale_factor=2, figure = figure)
for i in range(len(x_OB)):
	mlab.text3d(x_OB[i], y_OB[i], z_OB[i], OB_Names[i], scale=10, figure = figure)
#figure.scene.disable_render = False # Super duper trick
mlab.show()
