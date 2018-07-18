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
rel_path = "Data/A-OB-vrad-500pc-mean_AG_EBminR_v_xyz.fits"
abs_file_path = os.path.join(script_dir, rel_path)
rel_path2 = "Data/OB_positions.dat"
abs_file_path2 = os.path.join(script_dir, rel_path2)

# OB ass positions
ob = np.genfromtxt(abs_file_path2, names = True, dtype=None)
x_ob, y_ob, z_ob = ob['xg'], ob['yg'], ob['zg']
names_ob = ob['OB_Names']

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
#duration = 9
mlab.clf()
src = mlab.pipeline.scalar_field(b)
m = mlab.pipeline.iso_surface(src,
                              #contours = [0.2, 0.3, 0.5],
                              contours = [0.05,0.1,0.2,0.3],
                              opacity = 0.3, colormap = 'Blues', figure = fig_myv)
mlab.axes(#ranges = [-1300, 1300, -1300, 1300, -3500, 3500],
          ranges = [-454, 500, -497, 500, -273, 350],
          color = (0,0,0), nb_labels = 5, xlabel = 'X(pc)', ylabel = 'Y(pc)', zlabel = 'Z(pc)', figure = fig_myv)
mlab.points3d(x_ob_map, y_ob_map, z_ob_map, np.full(len(x_ob_map), 20),  scale_factor = 0.5)
for i in range(len(x_ob)):
    mlab.text3d(x_ob_map[i], y_ob_map[i], z_ob_map[i], str(names_ob[i]), scale = 5)
fig_myv.scene.disable_render = False
mlab.view(azimuth=0, elevation= 0,figure=fig_myv)
mlab.show()
