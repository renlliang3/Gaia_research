import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from sklearn.neighbors import KernelDensity
from mayavi import mlab
from  numpy.lib.recfunctions import append_fields
import astropy.coordinates as asc
import numpy as np
import random
import os

"""
# OB ass positions
ob = np.genfromtxt("Data/OB_positions.dat", names = True, dtype=None)
x_OB, y_OB, z_OB = ob['xg'], ob['yg'], ob['zg']
OB_Names = ob['OB_Names']
"""

readresults2 = Table.read('Data/OBRegions/A-OB-vrad-maxvtan40-mean_AG_EBminR/ScoCencontour7.5e10-8.dat',format='ascii')
results2 = np.array(readresults2)

x_OD, y_OD, z_OD, vx_OD, vy_OD, vz_OD = results2['xg'], results2['yg'], results2['zg'], results2['vx']+14.0, results2['vy']+12.24, results2['vz']+7.25

readresults = Table.read('Data/A-OB-vrad-500pc-maxvtan40-10vb-20-mean_AG_EBminR_v_xyz.fits',format='fits')
results = np.array(readresults)

#readresults2 = Table.read("Data/Updated-A-Query-Radial-Velocities.fits",format='fits')
#results2 = np.array(readresults2)

x, y, z, vx, vy, vz = results['xg'], results['yg'], results['zg'], results['vx'], results['vy'], results['vz']

xmin=-150
xmax=200
ymin=-275
ymax=75
zmin=-175
zmax=175

x_cube = x[np.where((x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)&(z>=zmin)&(z<=zmax))]
y_cube = y[np.where((x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)&(z>=zmin)&(z<=zmax))]
z_cube = z[np.where((x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)&(z>=zmin)&(z<=zmax))]

xyz = np.vstack([x_cube,y_cube,z_cube])

d = xyz.shape[0]
n = xyz.shape[1]
#bw = (n*(d+2)/4.)**(-1./(d+4))
bw = 20

kde_fit = KernelDensity(bandwidth=bw, kernel='gaussian').fit(xyz.T)

boxes=100j
X, Y, Z = np.mgrid[xmin:xmax:boxes, ymin:ymax:boxes, zmin:zmax:boxes]
positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), X.shape)

binsize = (xmax-xmin)/(boxes.imag)
w = int(boxes.imag/2)

figure = mlab.figure('myfig')
figure.scene.disable_render = True # Super duper trick
#mlab.contour3d(X, Y, Z, dens, extent = [-650, 650, -650, 650, -650, 650], opacity=0.3, colormap = 'Blues', figure = figure)
mlab.contour3d(X, Y, Z, dens, contours=[7.5*10**-8,8*10**-8,10*10**-8,12*10**-8], opacity=0.3, vmin=np.min(dens), vmax=np.max(dens), colormap = 'Blues', figure = figure)
#mlab.contour3d(X_OD, Y_OD, Z_OD, dens_OD, contours=[1*10**-8,2*10**-8,5*10**-8,8*10**-8,11*10**-8], opacity=0.3, vmin=np.min(dens_OD), vmax=np.max(dens_OD), colormap = 'Reds', figure = figure)
mlab.quiver3d(x_OD, y_OD, z_OD, vx_OD, vy_OD, vz_OD)
#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)
#mlab.points3d(x_OB, y_OB, z_OB, np.full(len(x_OB),10), scale_factor=2, figure = figure)
#mlab.points3d(x_cube, y_cube, z_cube, np.full(len(x_cube),1), scale_factor = 1, color=(0,0,0), figure = figure)
mlab.points3d(x_OD, y_OD, z_OD, np.full(len(x_OD),1), scale_factor = 1, color=(1,1,1), figure = figure)
mlab.axes(nb_labels = 5, figure = figure)
#for i in range(len(x_OB)):
#	mlab.text3d(x_OB[i], y_OB[i], z_OB[i], OB_Names[i], scale=10, figure = figure)

figure.scene.disable_render = False # Super duper trick
mlab.show()

"""
t = 14*10**-8
densmask = np.array((dens>=t), dtype = bool)

source_id_arr = np.array([],dtype='i8')
x_arr = np.array([],dtype='f8')
y_arr = np.array([],dtype='f8')
z_arr = np.array([],dtype='f8')
vx_arr = np.array([],dtype='f8')
vy_arr = np.array([],dtype='f8')
vz_arr = np.array([],dtype='f8')

for i in range(len(x)):
    #if x[i] >= xmin and x[i] <= xmax and y[i] >= ymin and y[i] <= ymax and z[i] >= zmin and z[i] <= zmax:
        x_voxel = int(round(float((x[i]-25)/binsize))) + w
        y_voxel = int(round(float((y[i]+100)/binsize))) + w
        z_voxel = int(round(float((z[i])/binsize))) + w
        if x_voxel >= 0 and x_voxel < 2*w and y_voxel >= 0 and y_voxel < 2*w and z_voxel >= 0 and z_voxel < 2*w and densmask[x_voxel][y_voxel][z_voxel] == True:
            source_id_arr = np.append(source_id_arr,results['source_id'][i])
            x_arr = np.append(x_arr,x[i])
            y_arr = np.append(y_arr,y[i])
            z_arr = np.append(z_arr,z[i])
            vx_arr = np.append(vx_arr,vx[i])
            vy_arr = np.append(vy_arr,vy[i])
            vz_arr = np.append(vz_arr,vz[i])

source_id_arr = np.rec.array(source_id_arr, dtype=[('source_id', np.int64)])
new_arr = append_fields(source_id_arr, 'xg', x_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'yg', y_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'zg', z_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'vx', vx_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'vy', vy_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'vz', vz_arr, usemask=False, dtypes=[np.float64])

np.savetxt('ScoCencontour14e10-8.dat', new_arr, header='source_id xg yg zg vx vy vz', fmt = '%2li %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f')
"""
