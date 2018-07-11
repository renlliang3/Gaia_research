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

readresults = Table.read("Data/Updated-A-Query-Radial-Velocities_v_xyz.fits",format='fits')
results = np.array(readresults)

readresults2 = Table.read("Data/Updated-A-Query-Radial-Velocities.fits",format='fits')
results2 = np.array(readresults2)

x = results['xg']
y = results['yg']
z = results['zg']

vx = results['vx']
vy = results['vy']
vz = results['vz']

v_xyz = np.vstack([vx,vy,vz])

d = v_xyz.shape[0]
n = v_xyz.shape[1]
#bw = (n*(d+2)/4.)**(-1./(d+4))
bw = 3

kde_fit = KernelDensity(bandwidth=bw, kernel='gaussian').fit(v_xyz.T)

rangesize=100
boxes=100j
VX, VY, VZ = np.mgrid[-rangesize:rangesize:boxes, -rangesize:rangesize:boxes, -rangesize:rangesize:boxes]
positions = np.vstack([VX.ravel(), VY.ravel(), VZ.ravel()])
dens = np.reshape(np.exp(kde_fit.score_samples(positions.T)), VX.shape)

binsize = 2*rangesize/(boxes.imag)
w = int(rangesize/binsize)
t = 6*10**-5
densmask = np.array((dens>t), dtype = bool)

source_id_arr = np.array([],dtype='i8')
x_arr = np.array([],dtype='f8')
y_arr = np.array([],dtype='f8')
z_arr = np.array([],dtype='f8')
vx_arr = np.array([],dtype='f8')
vy_arr = np.array([],dtype='f8')
vz_arr = np.array([],dtype='f8')

for i in range(len(x)):
	vx_voxel = int(round(float(vx[i]/binsize)))+w
    	vy_voxel = int(round(float(vy[i]/binsize)))+w
    	vz_voxel = int(round(float(vz[i]/binsize)))+w
    	if vx_voxel >= 0 and vx_voxel < 2*w and vy_voxel >= 0 and vy_voxel < 2*w and vz_voxel >= 0 and vz_voxel < 2*w and densmask[vx_voxel][vy_voxel][vz_voxel] == True:
		source_id_arr = np.append(source_id_arr,results2['source_id'][i])
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

np.savetxt('Sourceswithincontour.dat', new_arr, header='source_id xg yg zg vx vy vz', fmt = '%2li %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f')
