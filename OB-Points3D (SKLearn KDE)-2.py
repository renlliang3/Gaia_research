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
rel_path = "Data/A-OB-500pc-no-extinction_xyz.fits"
#rel_path = "Data/pms_500pc_C1C2_nogiants_prlx_2MASS_noext.fits"
abs_file_path = os.path.join(script_dir, rel_path)
rel_path2 = "Data/OB_positions.dat"
abs_file_path2 = os.path.join(script_dir, rel_path2)

readresults = Table.read(abs_file_path,format='fits')
results = np.array(readresults)

# OB ass positions
ob = np.genfromtxt(abs_file_path2, names = True, dtype=None)
x_OB, y_OB, z_OB = ob['xg'], ob['yg'], ob['zg']
OB_Names = ob['OB_Names']

x = results['xg']
y = results['yg']
z = results['zg']

#vx = coordinates_galactic.velocity.d_x.value
#vy = coordinates_galactic.velocity.d_y.value
#vz = coordinates_galactic.velocity.d_z.value

#vx = np.array([vx[i] for i in sorted(random_indices)])
#vy = np.array([vy[i] for i in sorted(random_indices)])
#vz = np.array([vz[i] for i in sorted(random_indices)])

xyz = np.vstack([x,y,z])

d = xyz.shape[0]
n = xyz.shape[1]
#bw = (n*(d+2)/4.)**(-1./(d+4))
bw = 20

kde_fit = KernelDensity(bandwidth=bw, kernel='gaussian').fit(xyz.T)

crange=500

X, Y, Z = np.mgrid[-crange:crange:100j, -crange:crange:100j, -crange:crange:100j]
positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), X.shape)

#Visualize the density estimate as isosurfaces
figure = mlab.figure('myfig')
figure.scene.disable_render = True # Super duper trick
#mlab.contour3d(X, Y, Z, dens, extent = [-650, 650, -650, 650, -650, 650], opacity=0.3, colormap = 'Blues', figure = figure)
mlab.contour3d(X, Y, Z, dens, contours=[8*10**-9,1*10**-8,1.5*10**-8,2*10**-8,2.5*10**-8,3*10**-8], opacity=0.3, vmin=np.min(dens), vmax=np.max(dens), colormap = 'Blues', figure = figure)
#contours=[10*10**-9,1.5*10**-8,3*10**-8,7*10**-8,1*10**-7],
#mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes(nb_labels = 5, figure = figure)
#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)
mlab.points3d(x_OB, y_OB, z_OB, np.full(len(x_OB),10), scale_factor=2, figure = figure)
for i in range(len(x_OB)):
	mlab.text3d(x_OB[i], y_OB[i], z_OB[i], OB_Names[i], scale=10, figure = figure)
figure.scene.disable_render = False # Super duper trick
mlab.show()
