import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from mayavi import mlab
import astropy.coordinates as asc
import numpy as np
import random
import os

script_dir = os.path.dirname(__file__)
rel_path = "Data/OBvrad1000-Katz_v_xyz.csv"
#rel_path = "Data/pms_500pc_C1C2_nogiants_prlx_2MASS_noext.fits"
abs_file_path = os.path.join(script_dir, rel_path)
rel_path2 = "Data/OB_positions.dat"
abs_file_path2 = os.path.join(script_dir, rel_path2)

readresults = Table.read(abs_file_path,format='csv')
results = np.array(readresults)

# OB ass positions
ob = np.genfromtxt(abs_file_path2, names = True, dtype=None)
x_OB, y_OB, z_OB = ob['xg'], ob['yg'], ob['zg']
OB_Names = ob['OB_Names']

x = results['xg']
y = results['yg']
z = results['zg']

#vx = coordinates_ICRS.velocity.d_x.value
#vy = coordinates_ICRS.velocity.d_y.value
#vz = coordinates_ICRS.velocity.d_z.value

#vx = coordinates_galactic.velocity.d_x.value
#vy = coordinates_galactic.velocity.d_y.value
#vz = coordinates_galactic.velocity.d_z.value

#random_indices = random.sample(range(len(x)),5000)

#x = np.array([x[i] for i in sorted(random_indices)])
#y = np.array([y[i] for i in sorted(random_indices)])
#z = np.array([z[i] for i in sorted(random_indices)])

#vx = np.array([vx[i] for i in sorted(random_indices)])
#vy = np.array([vy[i] for i in sorted(random_indices)])
#vz = np.array([vz[i] for i in sorted(random_indices)])

xyz = np.vstack([x,y,z])

kde = stats.gaussian_kde(xyz, bw_method=0.08)

crange=1000

X, Y, Z = np.mgrid[-crange:crange:100j, -crange:crange:100j, -crange:crange:100j]
positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
density = kde(positions).reshape(X.shape)

#Visualize the density estimate as isosurfaces
figure = mlab.figure('myfig')
figure.scene.disable_render = True # Super duper trick
#mlab.contour3d(X, Y, Z, dens, extent = [-650, 650, -650, 650, -650, 650], opacity=0.3, colormap = 'Blues', figure = figure)
mlab.contour3d(X, Y, Z, density, opacity=0.3, colormap = 'Blues', figure = figure)
#mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes(nb_labels = 5, figure = figure)
#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)
mlab.points3d(x_OB, y_OB, z_OB, np.full(len(x_OB),10), scale_factor=2, figure = figure)
for i in range(len(x_OB)):
	mlab.text3d(x_OB[i], y_OB[i], z_OB[i], OB_Names[i], scale=10, figure = figure)
figure.scene.disable_render = False # Super duper trick
mlab.show()
