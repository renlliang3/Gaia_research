import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf
from scipy.special import expit
import sys
import os
import imageio.core
from astropy.io import fits

#Config

#dataDir = '/data2/KevinJardine/tgas/'

#Code
binSize = 5
w = 650//binSize
z_height = 350//binSize


### Create map with the bright and young stars

sigma = 3
spread = 80
fp_hdu = fits.open('tgasHip_2M_hot.fits')
fp = fp_hdu[1].data



print('File is open')
xgal, ygal, zgal = fp['xg'], fp['yg'], fp['zg']

ra, dec = fp['RA'], fp['DEC'] 
plx, e_plx = fp['parallax'], fp['parallax_error']
pmra, e_pmra = fp['pmRA'], fp['pmRA_error']
pmdec, e_pmdec = fp['pmDE'], fp['pmDE_error']

gmag, jmag, jmag_error = fp['Gmag'], fp['Jmag'], fp['Jmag_error']
gflux, gflux_error = fp['flux_gmag'], fp['flux_gmag_error']
       

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

fig_myv = mlab.figure(size=(750,700), bgcolor=(1, 1, 1), fgcolor = (0,0,0))
duration = 9

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

#mlab.show()
animation = mpy.VideoClip(make_frame, duration=duration)
animation.write_gif("sinc2.gif", fps=15)
