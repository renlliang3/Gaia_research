import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from sklearn.neighbors import KernelDensity
from mayavi import mlab
from  numpy.lib.recfunctions import append_fields
import matplotlib.pyplot as plt
import astropy.coordinates as asc
import numpy as np
import random
import os


readresults_OD1 = Table.read('Data/OBRegions/A-OB-vrad-maxvtan40-mean_AG_EBminR/NGC3532contour3.5e10-7-Cluster1.dat',format='ascii')
results_OD1 = np.array(readresults_OD1)

x_OD1, y_OD1, z_OD1, vx_OD1, vy_OD1, vz_OD1 = results_OD1['xg'], results_OD1['yg'], results_OD1['zg'], results_OD1['vx']+14.0, results_OD1['vy']+12.24, results_OD1['vz']+7.25

readresults_OD2 = Table.read('Data/OBRegions/A-OB-vrad-maxvtan40-mean_AG_EBminR/NGC3532contour3.5e10-7-Cluster2.dat',format='ascii')
results_OD2 = np.array(readresults_OD2)

x_OD2, y_OD2, z_OD2, vx_OD2, vy_OD2, vz_OD2 = results_OD2['xg'], results_OD2['yg'], results_OD2['zg'], results_OD2['vx']+14.0, results_OD2['vy']+12.24, results_OD2['vz']+7.25

"""
readresults_OD3 = Table.read('Data/OBRegions/A-OB-vrad-maxvtan40-mean_AG_EBminR/ScoCencontour7.5e10-8-UVW-0.6e10-4-Cluster3.dat',format='ascii')
results_OD3 = np.array(readresults_OD3)

x_OD3, y_OD3, z_OD3, vx_OD3, vy_OD3, vz_OD3 = results_OD3['xg'], results_OD3['yg'], results_OD3['zg'], results_OD3['vx']+14.0, results_OD3['vy']+12.24, results_OD3['vz']+7.25

readresults_OD4 = Table.read('Data/OBRegions/A-OB-vrad-maxvtan40-mean_AG_EBminR/ScoCencontour7.5e10-8-UVW-0.6e10-4-Cluster4.dat',format='ascii')
results_OD4 = np.array(readresults_OD4)

x_OD4, y_OD4, z_OD4, vx_OD4, vy_OD4, vz_OD4 = results_OD4['xg'], results_OD4['yg'], results_OD4['zg'], results_OD4['vx']+14.0, results_OD4['vy']+12.24, results_OD4['vz']+7.25
"""

readresults = Table.read('Data/OBRegions/A-OB-vrad-maxvtan40-mean_AG_EBminR/NGC3532contour3.5e10-7.dat',format='ascii')
results = np.array(readresults)

source_id, x, y, z, vx, vy, vz = results['source_id'], results['xg'], results['yg'], results['zg'], results['vx']+14.0, results['vy']+12.24, results['vz']+7.25

v_xyz = np.vstack([vx,vy,vz])

bw = 3

kde_fit = KernelDensity(bandwidth=bw, kernel='gaussian').fit(v_xyz.T)

rangesize=100
boxes=100j

VX, VY, VZ = np.mgrid[-rangesize:rangesize:boxes, -rangesize:rangesize:boxes, -rangesize:rangesize:boxes]
positions = np.vstack([VX.ravel(), VY.ravel(), VZ.ravel()])
dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), VX.shape)

binsize = 2*rangesize/(boxes.imag)
w = int(boxes.imag/2)

#Visualize the density estimate as isosurfaces
figure = mlab.figure('myfig')
figure.scene.disable_render = True # Super duper trick
#mlab.contour3d(VX, VY, VZ, dens, contours=[0.25*10**-5,0.5*10**-5,1*10**-5,2*10**-5,2.5*10**-5,3*10**-5], opacity=0.3, colormap = 'Blues', figure = figure)
mlab.contour3d(VX, VY, VZ, dens, contours=[0.1*10**-4,0.5*10**-4,1*10**-4,1.5*10**-4,2*10**-4,2.5*10**-4], opacity=0.3, vmin=np.min(dens), vmax=np.max(dens), colormap = 'Reds', figure = figure)
mlab.axes(xlabel='U (km/s)', ylabel='V (km/s)', zlabel='W (km/s)', nb_labels = 3, figure = figure)
#mlab.quiver3d(x, y, z, vx, vy, vz)

#mlab.axes(extent = [-650, 650, -650, 650, -650, 650], ranges = [-650, 650, -650, 650, -650, 650], nb_labels = 7, figure = figure)
#mlab.points3d(vx, vy, vz, np.full(len(vx),1), scale_factor = 1, color=(1,1,1), figure = figure)
mlab.points3d(vx_OD1, vy_OD1, vz_OD1, np.full(len(x_OD1),1), scale_factor = 1, color=(0,1,0), figure = figure)
mlab.points3d(vx_OD2, vy_OD2, vz_OD2, np.full(len(x_OD2),1), scale_factor = 1, color=(0,0,1), figure = figure)
#mlab.points3d(vx_OD3, vy_OD3, vz_OD3, np.full(len(x_OD3),1), scale_factor = 1, color=(0,0,1), figure = figure)
#mlab.points3d(vx_OD4, vy_OD4, vz_OD4, np.full(len(x_OD4),1), scale_factor = 1, color=(1,1,1), figure = figure)
figure.scene.disable_render = False # Super duper trick
mlab.show()

"""
fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=True)

# We can set the number of bins with the `bins` kwarg
axs[0].hist(vx, bins=10)
axs[1].hist(vy, bins=10)
axs[2].hist(vz, bins=10)

xlabel = ['U (km/s)', 'V (km/s)', 'W (km/s)']
ticks = [range(-40,100,20), range(-50,50,25), range(-20,30,10)]
i=0
for ax in axs:
    ax.set(xlabel=xlabel[i], ylabel='Counts', xticks = ticks[i])
    i+=1

plt.savefig('1D-Histograms-UVW-NGC3532contour3.5e10-7.png')


t = 2.5*10**-4
densmask = np.array((dens>=t), dtype = bool)

source_id_arr = np.array([],dtype='i8')
x_arr = np.array([],dtype='f8')
y_arr = np.array([],dtype='f8')
z_arr = np.array([],dtype='f8')
vx_arr = np.array([],dtype='f8')
vy_arr = np.array([],dtype='f8')
vz_arr = np.array([],dtype='f8')

for i in range(len(x)):
    vx_voxel = int(round(float((vx[i])/binsize))) + w
    vy_voxel = int(round(float((vy[i])/binsize))) + w
    vz_voxel = int(round(float((vz[i])/binsize))) + w
    if vx_voxel >= 0 and vx_voxel < 2*w and vy_voxel >= 0 and vy_voxel < 2*w and vz_voxel >= 0 and vz_voxel < 2*w and densmask[vx_voxel][vy_voxel][vz_voxel] == True:
        source_id_arr = np.append(source_id_arr,source_id[i])
        x_arr = np.append(x_arr,x[i])
        y_arr = np.append(y_arr,y[i])
        z_arr = np.append(z_arr,z[i])
        vx_arr = np.append(vx_arr,vx[i]-14.0)
        vy_arr = np.append(vy_arr,vy[i]-12.24)
        vz_arr = np.append(vz_arr,vz[i]-7.25)

source_id_arr = np.rec.array(source_id_arr, dtype=[('source_id', np.int64)])
new_arr = append_fields(source_id_arr, 'xg', x_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'yg', y_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'zg', z_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'vx', vx_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'vy', vy_arr, usemask=False, dtypes=[np.float64])
new_arr = append_fields(new_arr, 'vz', vz_arr, usemask=False, dtypes=[np.float64])

np.savetxt('Orioncontour2.5e10-7-UVW-2.5e10-4.dat', new_arr, header='source_id xg yg zg vx vy vz', fmt = '%2li %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f')
"""
