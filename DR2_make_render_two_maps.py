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
rel_path = "Data/pms_500pc_xyz.fits"
abs_file_path = os.path.join(script_dir, rel_path)
rel_path2 = "Data/OB-Katz_xyz.fits"
abs_file_path2 = os.path.join(script_dir, rel_path2)
rel_path3 = "Data/OB_positions.dat"
abs_file_path3 = os.path.join(script_dir, rel_path3)

# OB ass positions
#fits_OB = fits.open(abs_file_path3)
#fits_OB = fits_OB[1].data
#x_ob, y_ob, z_ob = fits_OB['xg'], fits_OB['yg'], fits_OB['zg']
#names_ob = ['GammaVel', 'NGC2547', 'Trumpler10', 'NGC2516', 'VelaOB2', 'OrionOB1' , 'HIP22931' , 'HIP27345', 'HIP33175', 'NGC2232' , 'Bellatrix', 'Rigel' , 'Alnilam' , 'US', 'UCL', 'LCC' , 'IC2602' ,' aCar' , 'IC2391', 'HIP45189', 'NGC2451A' , 'Col135',  'HIP28944']

ob = np.genfromtxt(abs_file_path3, names = True, dtype=None)
x_ob, y_ob, z_ob = ob['xg'], ob['yg'], ob['zg']
names_ob = ob['OB_Names']
#Code
binSize = 3
w = 500//binSize
z_height = 350//binSize

### Create map with the bright and young stars

sigma_pms = 2
spread_pms = 30
fp_hdu_pms = fits.open(abs_file_path)
fp_pms = fp_hdu_pms[1].data

print('File is open')
xgal_pms, ygal_pms, zgal_pms = fp_pms['xg'], fp_pms['yg'], fp_pms['zg']

sigma_OB = 5
spread_OB = 50
fp_hdu_OB = fits.open(abs_file_path2)
fp_OB = fp_hdu_OB[1].data

xgal_OB, ygal_OB, zgal_OB = fp_OB['xg'], fp_OB['yg'], fp_OB['zg']
#vx, vy, vz = fp['vx'], fp['vy'], fp['vz']

#ra, dec = fp['RA'], fp['DEC'] 
#plx, e_plx = fp['parallax'], fp['parallax_error']
#pmra, e_pmra = fp['pmRA'], fp['pmRA_error']
#pmdec, e_pmdec = fp['pmDE'], fp['pmDE_error']

#gmag, jmag, jmag_error = fp['Gmag'], fp['Jmag'], fp['Jmag_error']
#gflux, gflux_error = fp['flux_gmag'], fp['flux_gmag_error']
       

a_pms = np.zeros((w*2,w*2,2*z_height))
print('Making Map')
N_pms = len(xgal_pms)
for  i in range(N_pms):
    x_pms = int(round(float(xgal_pms[i])/binSize))+w
    y_pms = int(round(float(ygal_pms[i])/binSize))+w
    z_pms = int(round(float(zgal_pms[i])/binSize))+z_height
    if x_pms >= 0 and x_pms < 2*w and y_pms >= 0 and y_pms < 2*w and z_pms >= 0 and z_pms < 2*z_height:
        a_pms[x_pms][y_pms][z_pms] += 1

a_OB = np.zeros((w*2,w*2,2*z_height))
N_OB = len(xgal_OB)
for  i in range(N_OB):
    x_OB = int(round(float(xgal_OB[i])/binSize))+w
    y_OB = int(round(float(ygal_OB[i])/binSize))+w
    z_OB = int(round(float(zgal_OB[i])/binSize))+z_height
    if x_OB >= 0 and x_OB < 2*w and y_OB >= 0 and y_OB < 2*w and z_OB >= 0 and z_OB < 2*z_height:
        a_OB[x_OB][y_OB][z_OB] += 1
	
#x_ob_map, y_ob_map, z_ob_map = x_ob, y_ob, z_ob
x_ob_map, y_ob_map, z_ob_map = np.full(len(x_ob),0), np.full(len(y_ob),0), np.full(len(z_ob),0)

for o in range(len(x_ob)):
    x_ob_map[o] = int(round(float(x_ob[o])/binSize))+w
    y_ob_map[o] = int(round(float(y_ob[o])/binSize))+w
    z_ob_map[o] = int(round(float(z_ob[o])/binSize))+z_height

print('Smoothing')
gaussian_pms = gf(a_pms, sigma=sigma_pms, truncate=3)
b_pms = 2*(expit(spread_pms*gaussian_pms)-0.5)

gaussian_OB = gf(a_OB, sigma=sigma_OB, truncate=3)
b_OB = 2*(expit(spread_OB*gaussian_OB)-0.5)

from mayavi import mlab
import  moviepy.editor as mpy

fig_myv = mlab.figure(size=(750,700), bgcolor=(1, 1, 1), fgcolor = (0,0,0))
fig_myv.scene.disable_render = True
#duration = 9
mlab.clf()
src_pms = mlab.pipeline.scalar_field(b_pms)
src_OB = mlab.pipeline.scalar_field(b_OB)
m_pms = mlab.pipeline.iso_surface(src_pms,
                              #contours = [0.2, 0.3, 0.5],
                              contours = [0.95, 0.97, 0.99],
                              opacity = 0.3, colormap = 'Reds', figure = fig_myv)
m_OB = mlab.pipeline.iso_surface(src_OB,
                              #contours = [0.2, 0.3, 0.5],
                              contours = [0.1,0.2,0.3],
                              opacity = 0.3, colormap = 'Blues', figure = fig_myv)
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
