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
rel_path = "Data/OBvrad650-Katz.csv"
abs_file_path = os.path.join(script_dir, rel_path)

readresults = Table.read(abs_file_path,format='csv')
results = np.array(readresults)
distances = 1000/results['parallax']

#Convert coordinates to galactic and then to cartesian. 
#coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, pm_ra_cosdec=results['pmra']*u.mas/u.yr, pm_dec=results['pmdec']*u.mas/u.yr, radial_velocity=results['radial_velocity']*u.km/u.s, frame='icrs', obstime='J2015.5')
coordinates_ICRS = asc.SkyCoord(ra=results['ra']*u.degree, dec=results['dec']*u.degree, distance=distances*u.pc, frame='icrs', obstime='J2015.5')

#coordinates_cartesian_2 = np.column_stack((coordinates_ICRS.cartesian.x.value, coordinates_ICRS.cartesian.y.value, coordinates_ICRS.cartesian.z.value))
coordinates_galactic = coordinates_ICRS.galactic
coordinates_cartesian = np.column_stack((coordinates_galactic.cartesian.x.value, coordinates_galactic.cartesian.y.value, coordinates_galactic.cartesian.z.value))

gammavelorum_ICRS = asc.SkyCoord(ra='08h09m31.95013s', dec=(-47-(20/60)-(11.7108/3600))*u.degree, distance=336*u.pc, frame='icrs')
gammavelorum_galactic = gammavelorum_ICRS.galactic
NGC2547_ICRS = asc.SkyCoord(ra=360*(8/24 + 10.7/(24*60))*u.degree,dec=-(49+16/60)*u.degree, distance=410*u.pc, frame='icrs')
NGC2547_galactic = NGC2547_ICRS.galactic
trumpler10_ICRS = asc.SkyCoord(ra=360*(8/24+47.8/(24*60))*u.degree, dec=-(42+29/60)*u.degree, distance=366*u.pc, frame='icrs')
trumpler10_galactic = trumpler10_ICRS.galactic
NGC2516_ICRS = asc.SkyCoord(ra=360*(7/24+58/(24*60)+20/(24*3600))*u.degree, dec=-(60+52/60)*u.degree, distance=373*u.pc, frame='icrs')
NGC2516_galactic = NGC2516_ICRS.galactic
OriOB1_galactic = asc.SkyCoord(l=206.91*u.degree, b=-17.59*u.degree, distance=400*u.pc, frame='galactic')

HIP22931_ICRS = asc.SkyCoord(ra=360*(4/24+55/(24*60)+14/(24*3600))*u.degree, dec=-(5+9.8/60)*u.degree, distance=180*u.pc, frame='icrs')
HIP22931_galactic = HIP22931_ICRS.galactic
HIP27345_ICRS = asc.SkyCoord(ra=360*(5/24+48/(24*60)+16/(24*3600))*u.degree, dec=-(4+0.7/60)*u.degree, distance=202*u.pc, frame='icrs')
HIP27345_galactic = HIP27345_ICRS.galactic
HIP33175_ICRS = asc.SkyCoord(ra=360*(6/24+57/(24*60)+14/(24*3600))*u.degree, dec=-(5+55.2/60)*u.degree, distance=195*u.pc, frame='icrs')
HIP33175_galactic = HIP33175_ICRS.galactic
NGC2232_ICRS = asc.SkyCoord(ra=360*(6/24+27/(24*60)+15/(24*3600))*u.degree, dec=-(4+45/60+30/3600)*u.degree, distance=325*u.pc, frame='icrs')
NGC2232_galactic = NGC2232_ICRS.galactic

Bellatrix_ICRS = asc.SkyCoord(ra=360*(5/24+25/(24*60)+7.86325/(24*3600))*u.degree, dec=+(6+20/60+58.9318/3600)*u.degree, distance=77*u.pc, frame='icrs')
Bellatrix_galactic = Bellatrix_ICRS.galactic
Rigel_ICRS = asc.SkyCoord(ra=360*(5/24+14/(24*60)+32.27/(24*3600))*u.degree, dec=-(8+12/60+5.8981/3600)*u.degree, distance=260*u.pc, frame='icrs')
Rigel_galactic = Rigel_ICRS.galactic
Alnilam_ICRS = asc.SkyCoord(ra=360*(5/24+36/(24*60)+12.8/(24*3600))*u.degree, dec=-(1+12/60+6.9/3600)*u.degree, distance=260*u.pc, frame='icrs')
Alnilam_galactic = Alnilam_ICRS.galactic

USco_ICRS = asc.SkyCoord(ra=360*(16/24+12/(24*60)+0/(24*3600))*u.degree, dec=-(23+24/60+0/3600)*u.degree, distance=145*u.pc, frame='icrs')
USco_galactic = USco_ICRS.galactic
UCL_ICRS = asc.SkyCoord(ra=360*(15/24+24/(24*60)+0/(24*3600))*u.degree, dec=-(41+54/60+0/3600)*u.degree, distance=140*u.pc, frame='icrs')
UCL_galactic = UCL_ICRS.galactic
LLC_ICRS = asc.SkyCoord(ra=360*(12/24+19/(24*60)+0/(24*3600))*u.degree, dec=-(57+06/60+0/3600)*u.degree, distance=118*u.pc, frame='icrs')
LLC_galactic = LLC_ICRS.galactic
IC2602_ICRS = asc.SkyCoord(ra=360*(10/24+42/(24*60)+58/(24*3600))*u.degree, dec=-(64+24/60+0/3600)*u.degree, distance=147*u.pc, frame='icrs')
IC2602_galactic = IC2602_ICRS.galactic
aCar_ICRS = asc.SkyCoord(ra=360*(9/24+9/(24*60)+30/(24*3600))*u.degree, dec=-(59+07/60+42/3600)*u.degree, distance=150*u.pc, frame='icrs')
aCar_galactic = aCar_ICRS.galactic
IC2391_ICRS = asc.SkyCoord(ra=360*(8/24+40/(24*60)+32/(24*3600))*u.degree, dec=-(53+02/60+0/3600)*u.degree, distance=176*u.pc, frame='icrs')
IC2391_galactic = IC2391_ICRS.galactic
HIP45189_ICRS = asc.SkyCoord(ra=360*(9/24+13/(24*60)+47/(24*3600))*u.degree, dec=-(43+44/60+24/3600)*u.degree, distance=174*u.pc, frame='icrs')
HIP45189_galactic = HIP45189_ICRS.galactic
NGC2451A_ICRS = asc.SkyCoord(ra=360*(7/24+43/(24*60)+12/(24*3600))*u.degree, dec=-(38+24/60+0/3600)*u.degree, distance=188*u.pc, frame='icrs')
NGC2451A_galactic = NGC2451A_ICRS.galactic
Col135_ICRS = asc.SkyCoord(ra=360*(7/24+17/(24*60)+17/(24*3600))*u.degree, dec=-(36+49/60+0/3600)*u.degree, distance=258*u.pc, frame='icrs')
Col135_galactic = Col135_ICRS.galactic
HIP28944_ICRS = asc.SkyCoord(ra=360*(6/24+9/(24*60)+36/(24*3600))*u.degree, dec=-(22+9.1/60+0/3600)*u.degree, distance=272*u.pc, frame='icrs')
HIP28944_galactic = HIP28944_ICRS.galactic

def addOB(x,y,z,*OBassociation):
	for item in OBassociation:
		x_item = item.cartesian.x.value
		y_item = item.cartesian.y.value
		z_item = item.cartesian.z.value
		x = np.append(x, x_item)
		y = np.append(y, y_item)
		z = np.append(z, z_item)
	return x,y,z

x_OB, y_OB, z_OB = np.array([]), np.array([]), np.array([])
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, gammavelorum_galactic, NGC2547_galactic, trumpler10_galactic, NGC2516_galactic, OriOB1_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, HIP22931_galactic, HIP27345_galactic, HIP33175_galactic, NGC2232_galactic, Bellatrix_galactic, Rigel_galactic, Alnilam_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, USco_galactic, UCL_galactic, LLC_galactic, IC2602_galactic, aCar_galactic, IC2391_galactic, HIP45189_galactic, NGC2451A_galactic, Col135_galactic, HIP28944_galactic)


#x_gammavelorum = gammavelorum_galactic.cartesian.x.value
#y_gammavelorum = gammavelorum_galactic.cartesian.y.value
#z_gammavelorum = gammavelorum_galactic.cartesian.z.value

#x_NGC2547 = NGC2547_galactic.cartesian.x.value
#y_NGC2547 = NGC2547_galactic.cartesian.y.value
#z_NGC2547 = NGC2547_galactic.cartesian.z.value

#x_OB = np.append(x_gammavelorum, x_NGC2547)
#y_OB = np.append(y_gammavelorum, y_NGC2547)
#z_OB = np.append(z_gammavelorum, z_NGC2547)

#Check if their galactic coordinates matches ours 
#realgalactic = np.column_stack((results['l'],results['b'],distances))
#convertedgalactic = np.column_stack((coordinates_galactic.l.degree,coordinates_galactic.b.degree,coordinates_galactic.distance.pc))
#differencegalactic = convertedgalactic-realgalactic

#Check for differences in the cartesian coordinates between using the converted and the real galactic coordinates
#coordinates_cartesian_galactic_real = asc.spherical_to_cartesian(distances*u.pc, results['b']*u.degree, results['l']*u.degree)
#coordinates_cartesian_galactic_real = np.column_stack((coordinates_cartesian_galactic_real[0].value, coordinates_cartesian_galactic_real[1].value, coordinates_cartesian_galactic_real[2].value))
#differencecartesian = coordinates_cartesian_galactic_real - coordinates_cartesian_galactic_converted

x = coordinates_cartesian[:,0]#.filled(0)
y = coordinates_cartesian[:,1]#.filled(0)
z = coordinates_cartesian[:,2]#.filled(0)

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
mlab.contour3d(X, Y, Z, density, opacity=0.5)
#mlab.quiver3d(x, y, z, vx, vy, vz)
mlab.axes()
mlab.points3d(x_OB, y_OB, z_OB, np.full(len(x_OB),20), scale_factor=1)
OB_Names = ['gammavelorum', 'NGC2547', 'trumpler10', 'NGC2516', 'OriOB1', 'HIP22931', 'HIP27345', 'HIP33175', 'NGC2232', 'Bellatrix', 'Rigel', 'Alnilam', 'USco', 'UCL', 'LLC', 'IC2602', 'aCar', 'IC2391', 'HIP45189', 'NGC2451A', 'Col135', 'HIP28944']
for i in range(len(x_OB)):
	mlab.text3d(x_OB[i], y_OB[i], z_OB[i], OB_Names[i], scale=10)
figure.scene.disable_render = False # Super duper trick
mlab.show()
