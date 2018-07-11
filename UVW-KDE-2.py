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

script_dir = os.path.dirname(__file__)
rel_path = "Data/OBvrad650-Katz_v_xyz.fits"
abs_file_path = os.path.join(script_dir, rel_path)

readresults = Table.read(abs_file_path,format='fits')
results = np.array(readresults)

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

VX, VY, VZ = np.mgrid[-100:100:100j, -100:100:100j, -50:50:100j]
positions = np.vstack([VX.ravel(), VY.ravel(), VZ.ravel()])
dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), VX.shape)

#Visualize the density estimate as isosurfaces
figure = mlab.figure('myfig')
figure.scene.disable_render = True # Super duper trick
mlab.contour3d(VX, VY, VZ, dens, opacity=0.3, colormap = 'Blues', figure = figure)
#mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes(xlabel='vx (km/s)', ylabel='vy (km/s)', zlabel='vz (km/s)', nb_labels = 5, figure = figure)
#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)

figure.scene.disable_render = False # Super duper trick
mlab.show()
