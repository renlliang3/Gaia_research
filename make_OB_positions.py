import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf
from scipy.special import expit
import sys
import os
import csv
import imageio.core
from astropy.io import fits

#Config

script_dir = os.path.dirname(__file__)
rel_path = "Data/OB-Katz_xyz.fits"
abs_file_path = os.path.join(script_dir, rel_path)

#Code
binSize = 5
w = 650//binSize
z_height = 350//binSize


### Create map with the bright and young stars

sigma = 3
spread = 80
fp_hdu = fits.open(abs_file_path)
fp = fp_hdu[1].data

print('File is open')
xgal, ygal, zgal = fp['xg'], fp['yg'], fp['zg']
vx, vy, vz = fp['vx'], fp['vy'], fp['vz']

#ra, dec = fp['RA'], fp['DEC'] 
#plx, e_plx = fp['parallax'], fp['parallax_error']
#pmra, e_pmra = fp['pmRA'], fp['pmRA_error']
#pmdec, e_pmdec = fp['pmDE'], fp['pmDE_error']

#gmag, jmag, jmag_error = fp['Gmag'], fp['Jmag'], fp['Jmag_error']
#gflux, gflux_error = fp['flux_gmag'], fp['flux_gmag_error']
       

a = np.zeros((w*2,w*2,2*z_height))
print('Making Map')
N = len(xgal)
for  i in range(N):
    x = int(round(float(xgal[i])/binSize))+w
    y = int(round(float(ygal[i])/binSize))+w
    z = int(round(float(zgal[i])/binSize))+z_height
    if x >= 0 and x < 2*w and y >= 0 and y < 2*w and z >= 0 and z < 2*z_height:
        a[x][y][z] += 1
	


print('Smoothing')
gaussian = gf(a, sigma=sigma, truncate=3)
b = 2*(expit(spread*gaussian)-0.5)

from mayavi import mlab
import  moviepy.editor as mpy
import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
import astropy.coordinates as asc

#gammavelorum_ICRS = asc.SkyCoord(ra='08h09m31.95013s', dec=(-47-(20/60)-(11.7108/3600))*u.degree, distance=336*u.pc, frame='icrs')
#gammavelorum_galactic = gammavelorum_ICRS.galactic
gammavelorum_galactic = asc.SkyCoord(l=262.8025*u.degree, b=-7.6858*u.degree, distance=336*u.pc, frame='galactic')
#NGC2547_ICRS = asc.SkyCoord(ra=360*(8/24 + 10.7/(24*60))*u.degree,dec=-(49+16/60)*u.degree, distance=410*u.pc, frame='icrs')
#NGC2547_galactic = NGC2547_ICRS.galactic
NGC2547_galactic = asc.SkyCoord(l=264.465*u.degree, b=-8.597*u.degree, distance=410*u.pc, frame='galactic')
#trumpler10_ICRS = asc.SkyCoord(ra=360*(8/24+47.8/(24*60))*u.degree, dec=-(42+29/60)*u.degree, distance=366*u.pc, frame='icrs')
#trumpler10_galactic = trumpler10_ICRS.galactic
trumpler10_galactic = asc.SkyCoord(l=262.791*u.degree, b=0.674*u.degree, distance=366*u.pc, frame='galactic')

#NGC2516_ICRS = asc.SkyCoord(ra=360*(7/24+58/(24*60)+20/(24*3600))*u.degree, dec=-(60+52/60)*u.degree, distance=373*u.pc, frame='icrs')
#NGC2516_galactic = NGC2516_ICRS.galactic
NGC2516_galactic = asc.SkyCoord(l=273.816*u.degree, b=-15.856*u.degree, distance=373*u.pc, frame='galactic')
velaOB2_galactic = asc.SkyCoord(l=263*u.degree, b=-7*u.degree, distance=410*u.pc, frame='galactic')
OriOB1_galactic = asc.SkyCoord(l=206.91*u.degree, b=-17.59*u.degree, distance=400*u.pc, frame='galactic')

#HIP22931_ICRS = asc.SkyCoord(ra=360*(4/24+55/(24*60)+14/(24*3600))*u.degree, dec=-(5+9.8/60)*u.degree, distance=180*u.pc, frame='icrs')
#HIP22931_galactic = HIP22931_ICRS.galactic
#HIP27345_ICRS = asc.SkyCoord(ra=360*(5/24+48/(24*60)+16/(24*3600))*u.degree, dec=-(4+0.7/60)*u.degree, distance=202*u.pc, frame='icrs')
#HIP27345_galactic = HIP27345_ICRS.galactic
#HIP33175_ICRS = asc.SkyCoord(ra=360*(6/24+57/(24*60)+14/(24*3600))*u.degree, dec=-(5+55.2/60)*u.degree, distance=195*u.pc, frame='icrs')
#HIP33175_galactic = HIP33175_ICRS.galactic
#NGC2232_ICRS = asc.SkyCoord(ra=360*(6/24+27/(24*60)+15/(24*3600))*u.degree, dec=-(4+45/60+30/3600)*u.degree, distance=325*u.pc, frame='icrs')
#NGC2232_galactic = NGC2232_ICRS.galactic
NGC2232_galactic = asc.SkyCoord(l=214.432*u.degree, b=-7.542*u.degree, distance=325*u.pc, frame='galactic')

#Bellatrix_ICRS = asc.SkyCoord(ra=360*(5/24+25/(24*60)+7.86325/(24*3600))*u.degree, dec=+(6+20/60+58.9318/3600)*u.degree, distance=77*u.pc, frame='icrs')
#Bellatrix_galactic = Bellatrix_ICRS.galactic
#Rigel_ICRS = asc.SkyCoord(ra=360*(5/24+14/(24*60)+32.27/(24*3600))*u.degree, dec=-(8+12/60+5.8981/3600)*u.degree, distance=260*u.pc, frame='icrs')
#Rigel_galactic = Rigel_ICRS.galactic
#Alnilam_ICRS = asc.SkyCoord(ra=360*(5/24+36/(24*60)+12.8/(24*3600))*u.degree, dec=-(1+12/60+6.9/3600)*u.degree, distance=260*u.pc, frame='icrs')
#Alnilam_galactic = Alnilam_ICRS.galactic

#USco_ICRS = asc.SkyCoord(ra=360*(16/24+12/(24*60)+0/(24*3600))*u.degree, dec=-(23+24/60+0/3600)*u.degree, distance=145*u.pc, frame='icrs')
#USco_galactic = USco_ICRS.galactic
USco_galactic = asc.SkyCoord(l=315.5055*u.degree, b=20.0228*u.degree, distance=145*u.pc, frame='galactic')
#UCL_ICRS = asc.SkyCoord(ra=360*(15/24+24/(24*60)+0/(24*3600))*u.degree, dec=-(41+54/60+0/3600)*u.degree, distance=140*u.pc, frame='icrs')
#UCL_galactic = UCL_ICRS.galactic
UCL_galactic = asc.SkyCoord(l=331.0211*u.degree, b=12.5023*u.degree, distance=140*u.pc, frame='galactic')

#LCC_ICRS = asc.SkyCoord(ra=360*(12/24+19/(24*60)+0/(24*3600))*u.degree, dec=-(57+6/60+0/3600)*u.degree, distance=118*u.pc, frame='icrs')
#LCC_galactic = LCC_ICRS.galactic
LCC_galactic = asc.SkyCoord(l=298.5171*u.degree, b=5.4934*u.degree, distance=118*u.pc, frame='galactic')

#IC2602_ICRS = asc.SkyCoord(ra=360*(10/24+42/(24*60)+58/(24*3600))*u.degree, dec=-(64+24/60+0/3600)*u.degree, distance=147*u.pc, frame='icrs')
#IC2602_galactic = IC2602_ICRS.galactic
IC2602_galactic = asc.SkyCoord(l=289.601*u.degree, b=-4.906*u.degree, distance=147*u.pc, frame='galactic')

#aCar_ICRS = asc.SkyCoord(ra=360*(9/24+9/(24*60)+30/(24*3600))*u.degree, dec=-(59+7/60+42/3600)*u.degree, distance=150*u.pc, frame='icrs')
#aCar_galactic = aCar_ICRS.galactic
aCar_galactic = asc.SkyCoord(l=277.6823*u.degree, b=-7.6209*u.degree, distance=150*u.pc, frame='galactic')

#IC2391_ICRS = asc.SkyCoord(ra=360*(8/24+40/(24*60)+32/(24*3600))*u.degree, dec=-(53+2/60+0/3600)*u.degree, distance=176*u.pc, frame='icrs')
#IC2391_galactic = IC2391_ICRS.galactic
IC2391_galactic = asc.SkyCoord(l=270.362*u.degree, b=-6.839*u.degree, distance=176*u.pc, frame='galactic')

#HIP45189_ICRS = asc.SkyCoord(ra=360*(9/24+13/(24*60)+47/(24*3600))*u.degree, dec=-(43+44/60+24/3600)*u.degree, distance=174*u.pc, frame='icrs')
#HIP45189_galactic = HIP45189_ICRS.galactic
#NGC2451A_ICRS = asc.SkyCoord(ra=360*(7/24+43/(24*60)+12/(24*3600))*u.degree, dec=-(38+24/60+0/3600)*u.degree, distance=188*u.pc, frame='icrs')
#NGC2451A_galactic = NGC2451A_ICRS.galactic
NGC2451A_galactic = asc.SkyCoord(l=252.575*u.degree, b=-7.299*u.degree, distance=188*u.pc, frame='galactic')

#Col135_ICRS = asc.SkyCoord(ra=360*(7/24+17/(24*60)+17/(24*3600))*u.degree, dec=-(36+49/60+0/3600)*u.degree, distance=258*u.pc, frame='icrs')
#Col135_galactic = Col135_ICRS.galactic
Col135_galactic = asc.SkyCoord(l=248.765*u.degree, b=-11.132*u.degree, distance=258*u.pc, frame='galactic')
#HIP28944_ICRS = asc.SkyCoord(ra=360*(6/24+9/(24*60)+36/(24*3600))*u.degree, dec=-(22+9.1/60+0/3600)*u.degree, distance=272*u.pc, frame='icrs')
#HIP28944_galactic = HIP28944_ICRS.galactic

Stock2_galactic = asc.SkyCoord(l=133.334*u.degree, b=-1.694*u.degree, distance=303*u.pc, frame='galactic')
IC348_galactic = asc.SkyCoord(l=160.490*u.degree, b=-17.802*u.degree, distance=385*u.pc, frame='galactic')
Stock10_galactic = asc.SkyCoord(l=171.62*u.degree, b=3.55*u.degree, distance=380*u.pc, frame='galactic')
NGC6405_galactic = asc.SkyCoord(l=356.580*u.degree, b=-0.777*u.degree, distance=487*u.pc, frame='galactic')
NGC7092_galactic = asc.SkyCoord(l=92.403*u.degree, b=-2.242*u.degree, distance=326*u.pc, frame='galactic')
IC4756_galactic = asc.SkyCoord(l=36.381*u.degree, b=5.242*u.degree, distance=484*u.pc, frame='galactic')
Stock1_galactic = asc.SkyCoord(l=60.247*u.degree, b=2.257*u.degree, distance=320*u.pc, frame='galactic')
Pleiades_galactic = asc.SkyCoord(l=166.571*u.degree, b=-23.521*u.degree, distance=136*u.pc, frame='galactic')
NGC6475_galactic = asc.SkyCoord(l=355.861*u.degree, b=-4.501*u.degree, distance=301*u.pc, frame='galactic')
AlphaPer_galactic = asc.SkyCoord(l=146.568*u.degree, b=-5.862*u.degree, distance=185*u.pc, frame='galactic')
PerOB2_galactic = asc.SkyCoord(l=159.2406*u.degree, b=-17.1388*u.degree, distance=300*u.pc, frame='galactic')

NGC2422_galactic = asc.SkyCoord(l=230.958*u.degree, b=3.130*u.degree, distance=490*u.pc, frame='galactic')
Alessi21_galactic = asc.SkyCoord(l=223.44*u.degree, b=0.00*u.degree, distance=500*u.pc, frame='galactic')
Haffner13_galactic = asc.SkyCoord(l=245.00*u.degree, b=-3.73*u.degree, distance=714*u.pc, frame='galactic')
NGC2527_galactic = asc.SkyCoord(l=246.087*u.degree, b=1.855*u.degree, distance=601*u.pc, frame='galactic')
IC2395_galactic = asc.SkyCoord(l=266.637*u.degree, b=-3.592*u.degree, distance=705*u.pc, frame='galactic')
NGC3532_galactic = asc.SkyCoord(l=289.571*u.degree, b=1.347*u.degree, distance=486*u.pc, frame='galactic')

Roslund6_galactic = asc.SkyCoord(l=78.06*u.degree, b=0.26*u.degree, distance=450*u.pc, frame='galactic')
Platais1_galactic = asc.SkyCoord(l=92.561*u.degree, b=-1.646*u.degree, distance=1268*u.pc, frame='galactic')
Ascc122_galactic = asc.SkyCoord(l=95.91*u.degree, b=-15.9*u.degree, distance=700*u.pc, frame='galactic')
Stock12_galactic = asc.SkyCoord(l=111.44*u.degree, b=-8.48*u.degree, distance=400*u.pc, frame='galactic')

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
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, gammavelorum_galactic, NGC2547_galactic, trumpler10_galactic, NGC2516_galactic, velaOB2_galactic, OriOB1_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, NGC2232_galactic, USco_galactic, UCL_galactic, LCC_galactic, IC2602_galactic, aCar_galactic, IC2391_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, NGC2451A_galactic, Col135_galactic, Stock2_galactic, IC348_galactic, Stock10_galactic, NGC6405_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, NGC7092_galactic, IC4756_galactic, Stock1_galactic, Pleiades_galactic, NGC6475_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, AlphaPer_galactic, PerOB2_galactic, NGC2422_galactic, Alessi21_galactic, Haffner13_galactic, NGC2527_galactic, IC2395_galactic, NGC3532_galactic, Roslund6_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, Platais1_galactic, Ascc122_galactic, Stock12_galactic)

x_ob_map, y_ob_map, z_ob_map = np.full(len(x_OB),0), np.full(len(y_OB),0), np.full(len(z_OB),0)

for o in range(len(x_OB)):
    x_ob_map[o] = int(round(float(x_OB[o])/binSize))+w
    y_ob_map[o] = int(round(float(y_OB[o])/binSize))+w
    z_ob_map[o] = int(round(float(z_OB[o])/binSize))+z_height

#fig_myv = mlab.figure('myfig')
fig_myv = mlab.figure(size=(750,700), bgcolor=(1, 1, 1), fgcolor = (0,0,0))
fig_myv.scene.disable_render = True

duration = 9
mlab.clf()
src = mlab.pipeline.scalar_field(b)
m = mlab.pipeline.iso_surface(src, contours = [0.6, 0.8, 1], opacity = 0.3, colormap = 'Blues', figure = fig_myv)
mlab.axes(ranges = [-650, 650, -650, 650, -350, 350], color = (0,0,0), nb_labels = 5, xlabel = 'X(pc)', ylabel = 'Y(pc)', zlabel = 'Z(pc)', figure = fig_myv)
mlab.view(azimuth=0, elevation= 0, distance= 500, figure=fig_myv)

mlab.points3d(x_ob_map, y_ob_map, z_ob_map, np.full(len(x_OB),10), scale_factor=0.5)
#OB_Names = ['gammavelorum', 'NGC2547', 'trumpler10', 'NGC2516', 'velaOB2', 'OriOB1', 'HIP22931', 'HIP27345', 'HIP33175', 'NGC2232', 'Bellatrix', 'Rigel', 'Alnilam', 'USco', 'UCL', 'LLC', 'IC2602', 'aCar', 'IC2391', 'HIP45189', 'NGC2451A', 'Col135', 'HIP28944', 'Stock2', 'IC348', 'Stock10', 'NGC6405', 'NGC7092', 'IC4756', 'Stock1', 'Pleiades', 'NGC6475', 'AlphaPer', 'PerOB2', 'NGC2422', 'Alessi21', 'Haffner13', 'NGC2527', 'IC2395', 'NGC3532', 'Roslund6', 'Platais1', 'Ascc122', 'Stock12']

OB_Names = ['GammaVel', 'NGC2547', 'Trumpler10', 'NGC2516', 'VelaOB2', 'OriOB1', 'NGC2232', 'USco', 'UCL', 'LCC', 'IC2602', 'aCar', 'IC2391', 'NGC2451A', 'Col135', 'Stock2', 'IC348', 'Stock10', 'NGC6405', 'NGC7092', 'IC4756', 'Stock1', 'Pleiades', 'NGC6475', 'AlphaPer', 'PerOB2', 'NGC2422', 'Alessi21', 'Haffner13', 'NGC2527', 'IC2395', 'NGC3532', 'Roslund6', 'Platais1', 'Ascc122', 'Stock12']

for i in range(len(x_ob_map)):
	mlab.text3d(x_ob_map[i], y_ob_map[i], z_ob_map[i], OB_Names[i], scale=5)

fig_myv.scene.disable_render = False
#mlab.show()
#animation = mpy.VideoClip(make_frame, duration=duration)
#animation.write_gif("sinc2.gif", fps=15)

towrite = np.column_stack((OB_Names,x_OB,y_OB,z_OB)).tolist()
towrite.insert(0,['OB_Names','xg','yg','zg'])

with open('OB_positions.csv', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
    wr.writerows(towrite)
