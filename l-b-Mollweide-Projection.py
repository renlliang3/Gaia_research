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
rel_path = "Data/OB-Katz.csv"
#rel_path = "Data/pms_500pc_C1C2_nogiants_prlx_2MASS_noext.fits"
abs_file_path = os.path.join(script_dir, rel_path)

readresults = Table.read(abs_file_path,format='csv')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='icrs', obstime='J2015.5')
#coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, frame='icrs', obstime='J2015.5')

#coordinates_cartesian_2 = np.column_stack((coordinates_ICRS.cartesian.x.value, coordinates_ICRS.cartesian.y.value, coordinates_ICRS.cartesian.z.value))
coordinates_galactic = coordinates_ICRS.galactic
#coordinates_galactic = asc.SkyCoord(l=results['l']*u.degree, b=results['b']*u.degree, distance=distances*u.pc, frame='galactic', obstime='J2015.5')
coordinates_cartesian = np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))


x = coordinates_cartesian[:,0]#.filled(0)
y = coordinates_cartesian[:,1]#.filled(0)
z = coordinates_cartesian[:,2]#.filled(0)

listindices = np.array([])

for i in range(len(x))
	if x[i]>-325 && y[i]>0:
		listindices=np.append(listindices,i)

		
list = [i for i in x if x

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

kde = stats.gaussian_kde(xyz, bw_method=0.1)

X, Y, Z = np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]
positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
density = kde(positions).reshape(X.shape)

#Visualize the density estimate as isosurfaces
figure = mlab.figure('myfig')
figure.scene.disable_render = True # Super duper trick
mlab.contour3d(X, Y, Z, density, opacity=0.3, colormap = 'Blues', figure = figure)
#mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)
mlab.points3d(x_OB, y_OB, z_OB, np.full(len(x_OB),10), scale_factor=1, figure = figure)
OB_Names = ['gammavelorum', 'NGC2547', 'trumpler10', 'NGC2516', 'velaOB2', 'OriOB1', 'HIP22931', 'HIP27345', 'HIP33175', 'NGC2232', 'Bellatrix', 'Rigel', 'Alnilam', 'USco', 'UCL', 'LLC', 'IC2602', 'aCar', 'IC2391', 'HIP45189', 'NGC2451A', 'Col135', 'HIP28944']
for i in range(len(x_OB)):
	mlab.text3d(x_OB[i], y_OB[i], z_OB[i], OB_Names[i], scale=5, figure = figure)
figure.scene.disable_render = False # Super duper trick
mlab.show()
