import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf
from scipy.special import expit
import sys
import os
import imageio.core
from astropy.io import fits

#Config

script_dir = os.path.dirname(__file__)
#rel_path = "Data/OBvrad650-Katz_xyz.fits"
rel_path = "Data/OB-Katz_xyz.fits"
abs_file_path = os.path.join(script_dir, rel_path)
rel_path2 = "Data/OB_positions.fits"
abs_file_path2 = os.path.join(script_dir, rel_path2)

# OB ass positions
fits_OB = fits.open(abs_file_path2)
fits_OB = fits_OB[1].data
#ob = np.genfromtxt('OB_positions', names = True)
x_ob, y_ob, z_ob = fits_OB['xg'], fits_OB['yg'], fits_OB['zg']

names_ob = ['GammaVel', 'NGC2547', 'Trumpler10', 'NGC2516', 'VelaOB2', 'OrionOB1' , 'HIP22931' , 'HIP27345', 'HIP33175',
                'NGC2232' , 'Bellatrix', 'Rigel' , 'Alnilam' , 'US', 'UCL', 'LCC' , 'IC2602' ,
                ' aCar' , 'IC2391', 'HIP45189', 'NGC2451A' , 'Col135',  'HIP28944']
#Code
binSize = 3
w = 500//binSize
z_height = 350//binSize


### Create map with the bright and young stars

sigma = 5
spread = 50
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
	
#x_ob_map, y_ob_map, z_ob_map = x_ob, y_ob, z_ob
x_ob_map, y_ob_map, z_ob_map = np.full(len(x_ob),0), np.full(len(y_ob),0), np.full(len(z_ob),0)

for o in range(len(x_ob)):
    x_ob_map[o] = int(round(float(x_ob[o])/binSize))+w
    y_ob_map[o] = int(round(float(y_ob[o])/binSize))+w
    z_ob_map[o] = int(round(float(z_ob[o])/binSize))+z_height

print('Smoothing')
gaussian = gf(a, sigma=sigma, truncate=3)
b = 2*(expit(spread*gaussian)-0.5)

from mayavi import mlab
import  moviepy.editor as mpy

fig_myv = mlab.figure(size=(750,700), bgcolor=(1, 1, 1), fgcolor = (0,0,0))
fig_myv.scene.disable_render = True
#duration = 9
mlab.clf()
src = mlab.pipeline.scalar_field(b)
m = mlab.pipeline.iso_surface(src,
                              #contours = [0.2, 0.3, 0.5],
                              contours = [0.1,0.2,.3],
                              opacity = 0.3, colormap = 'Reds', figure = fig_myv)
mlab.axes(#ranges = [-1300, 1300, -1300, 1300, -3500, 3500],
          ranges = [-500, 500, -500, 500, -350, 350],
          color = (0,0,0), nb_labels = 5, xlabel = 'X(pc)', ylabel = 'Y(pc)', zlabel = 'Z(pc)', figure = fig_myv)
mlab.points3d(x_ob_map, y_ob_map, z_ob_map, np.full(len(x_ob_map), 20),  scale_factor = 0.5)
for i in range(len(x_ob)):
    mlab.text3d(x_ob_map[i], y_ob_map[i], z_ob_map[i], str(names_ob[i]), scale = 5)
fig_myv.scene.disable_render = False
mlab.view(azimuth=0, elevation= 0,figure=fig_myv)


mlab.show()



#def make_frame(t):
#    mlab.clf()
#    src = mlab.pipeline.scalar_field(b)
#    m = mlab.pipeline.iso_surface(src, contours = [0.8, 0.9, 1], opacity = 0.3,
#                                  colormap = 'Blues', figure = fig_myv)
#    mlab.axes(ranges = [-650, 650, -650, 650, -350, 350], color = (0,0,0), nb_labels = 5,
#                      xlabel = 'X(pc)', ylabel = 'Y(pc)', zlabel = 'Z(pc)', figure = fig_myv)
#    mlab.view(azimuth=0, elevation= 0, distance= 500, figure=fig_myv)

    #m.scene.camera.elevation(-95)
        
#    m.scene.camera.elevation(-10*t)
#    m.scene.camera.orthogonalize_view_up()
    #m.scene.camera.compute_view_plane_normal()
#    return mlab.screenshot(antialiased = True)
#animation = mpy.VideoClip(make_frame, duration=duration)
#animation.write_gif("sinc2.gif", fps=15)
