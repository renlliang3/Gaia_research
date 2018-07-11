import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf
from scipy.special import expit
import sys
import os
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

gammavelorum_ICRS = asc.SkyCoord(ra='08h09m31.95013s', dec=(-47-(20/60)-(11.7108/3600))*u.degree, distance=336*u.pc, frame='icrs')
gammavelorum_galactic = gammavelorum_ICRS.galactic
NGC2547_ICRS = asc.SkyCoord(ra=360*(8/24 + 10.7/(24*60))*u.degree,dec=-(49+16/60)*u.degree, distance=410*u.pc, frame='icrs')
NGC2547_galactic = NGC2547_ICRS.galactic
trumpler10_ICRS = asc.SkyCoord(ra=360*(8/24+47.8/(24*60))*u.degree, dec=-(42+29/60)*u.degree, distance=366*u.pc, frame='icrs')
trumpler10_galactic = trumpler10_ICRS.galactic
NGC2516_ICRS = asc.SkyCoord(ra=360*(7/24+58/(24*60)+20/(24*3600))*u.degree, dec=-(60+52/60)*u.degree, distance=373*u.pc, frame='icrs')
NGC2516_galactic = NGC2516_ICRS.galactic
velaOB2_galactic = asc.SkyCoord(l=263*u.degree, b=-7*u.degree, distance=410*u.pc, frame='galactic')
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
LLC_ICRS = asc.SkyCoord(ra=360*(12/24+19/(24*60)+0/(24*3600))*u.degree, dec=-(57+6/60+0/3600)*u.degree, distance=118*u.pc, frame='icrs')
LLC_galactic = LLC_ICRS.galactic
IC2602_ICRS = asc.SkyCoord(ra=360*(10/24+42/(24*60)+58/(24*3600))*u.degree, dec=-(64+24/60+0/3600)*u.degree, distance=147*u.pc, frame='icrs')
IC2602_galactic = IC2602_ICRS.galactic
aCar_ICRS = asc.SkyCoord(ra=360*(9/24+9/(24*60)+30/(24*3600))*u.degree, dec=-(59+7/60+42/3600)*u.degree, distance=150*u.pc, frame='icrs')
aCar_galactic = aCar_ICRS.galactic
IC2391_ICRS = asc.SkyCoord(ra=360*(8/24+40/(24*60)+32/(24*3600))*u.degree, dec=-(53+2/60+0/3600)*u.degree, distance=176*u.pc, frame='icrs')
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
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, gammavelorum_galactic, NGC2547_galactic, trumpler10_galactic, NGC2516_galactic, velaOB2_galactic, OriOB1_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, HIP22931_galactic, HIP27345_galactic, HIP33175_galactic, NGC2232_galactic, Bellatrix_galactic, Rigel_galactic, Alnilam_galactic)
x_OB, y_OB, z_OB = addOB(x_OB, y_OB, z_OB, USco_galactic, UCL_galactic, LLC_galactic, IC2602_galactic, aCar_galactic, IC2391_galactic, HIP45189_galactic, NGC2451A_galactic, Col135_galactic, HIP28944_galactic)

#fig_myv = mlab.figure('myfig')
fig_myv = mlab.figure(size=(750,700), bgcolor=(1, 1, 1), fgcolor = (0,0,0))
fig_myv.scene.disable_render = True

duration = 9
mlab.clf()
src = mlab.pipeline.scalar_field(b)
m = mlab.pipeline.iso_surface(src, contours = [0.6, 0.8, 1], opacity = 0.3, colormap = 'Blues', figure = fig_myv)
mlab.axes(ranges = [-650, 650, -650, 650, -350, 350], color = (0,0,0), nb_labels = 5, xlabel = 'X(pc)', ylabel = 'Y(pc)', zlabel = 'Z(pc)', figure = fig_myv)
mlab.view(azimuth=0, elevation= 0, distance= 500, figure=fig_myv)

mlab.points3d(x_OB, y_OB, z_OB, np.full(len(x_OB),20), scale_factor=1)
OB_Names = ['gammavelorum', 'NGC2547', 'trumpler10', 'NGC2516', 'velaOB2', 'OriOB1', 'HIP22931', 'HIP27345', 'HIP33175', 'NGC2232', 'Bellatrix', 'Rigel', 'Alnilam', 'USco', 'UCL', 'LLC', 'IC2602', 'aCar', 'IC2391', 'HIP45189', 'NGC2451A', 'Col135', 'HIP28944']
for i in range(len(x_OB)):
	mlab.text3d(x_OB[i], y_OB[i], z_OB[i], OB_Names[i], scale=10)

def make_frame(t):
    mlab.clf()
    src = mlab.pipeline.scalar_field(b)
    m = mlab.pipeline.iso_surface(src, contours = [0.8, 0.9, 1], opacity = 0.3,
                                  colormap = 'Blues', figure = fig_myv)
    mlab.axes(ranges = [-650, 650, -650, 650, -350, 350], color = (0,0,0), nb_labels = 5,
                      xlabel = 'X(pc)', ylabel = 'Y(pc)', zlabel = 'Z(pc)', figure = fig_myv)
    mlab.view(azimuth=0, elevation= 0, distance= 500, figure=fig_myv)

    #m.scene.camera.elevation(-95)
        
    m.scene.camera.elevation(-10*t)
    m.scene.camera.orthogonalize_view_up()
    #m.scene.camera.compute_view_plane_normal()
    return mlab.screenshot(antialiased = True)

fig_myv.scene.disable_render = False
mlab.show()
#animation = mpy.VideoClip(make_frame, duration=duration)
#animation.write_gif("sinc2.gif", fps=15)
