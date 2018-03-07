import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
from matplotlib.colors import LogNorm
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import numpy as np

job = Gaia.launch_job_async("select gaia.source_id, gaia.hip, \
gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 as g_mag_abs_gaia, \
gaia.phot_g_mean_mag+5*log10(hip.plx)-10 as g_mag_abs_hip, \
hip.b_v \
from gaiadr1.tgas_source as gaia \
inner join public.hipparcos_newreduction as hip \
on gaia.hip = hip.hip \
where gaia.parallax/gaia.parallax_error >= 5 and \
hip.plx/hip.e_plx >= 5 and \
hip.e_b_v > 0.0 and hip.e_b_v <= 0.05 and \
2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux <= 0.05", dump_to_file=True)

print(job)

r = job.get_results()
x=r['b_v'].filled(0)
y_gaia=r['g_mag_abs_gaia'].filled(0)
y_hip=r['g_mag_abs_hip'].filled(0)

xy_gaia=np.vstack([x,y_gaia])
xy_hip=np.vstack([x,y_hip])

d_gaia=xy_gaia.shape[0]
d_hip=xy_hip.shape[0]
n_gaia=xy_gaia.shape[1]
n_hip=xy_hip.shape[1]

print(d_gaia)
print(n_gaia)
print(d_hip)
print(n_hip)

bw_gaia = (n_gaia * (d_gaia + 2) / 4.)**(-1. / (d_gaia + 4))
bw_hip = (n_hip * (d_hip + 2) / 4.)**(-1. / (d_hip + 4))

kde_gaia_fit = KernelDensity(bandwidth=bw_gaia, kernel='epanechnikov').fit(xy_gaia.T)
kde_hip_fit = KernelDensity(bandwidth=bw_hip, kernel='epanechnikov').fit(xy_hip.T) 



X_gaia, Y_gaia=np.mgrid[x.min():x.max():100j, y_gaia.min():y_gaia.max():100j]
positions_gaia=np.vstack([X_gaia.ravel(), Y_gaia.ravel()])
dens_gaia = np.reshape(np.exp(kde_gaia_fit.score_samples(positions_gaia.T)), X_gaia.shape)

X_hip, Y_hip=np.mgrid[x.min():x.max():100j, y_hip.min():y_hip.max():100j]
positions_hip=np.vstack([X_hip.ravel(), Y_hip.ravel()])
dens_hip = np.reshape(np.exp(kde_hip_fit.score_samples(positions_hip.T)), X_hip.shape)


plt.subplot(1, 2, 1)
plt.imshow(np.rot90(dens_hip), cmap=plt.cm.viridis, extent=[x.min(), x.max(), y_hip.min(), y_hip.max()], aspect='auto')
#plt.contour(X_hip, Y_hip, dens_hip, colors='k', linewidth=0.01, levels = [10,30,50])
plt.text(-0.1, 10.0, 'Hipparcos')
plt.xlim(-0.3,2.0)
#plt.ylim(-4.0,12.5)
plt.xlabel('(B-V)')
plt.ylabel(r'$M_G$')
plt.gca().invert_yaxis()

plt.subplot(1, 2, 2)
plt.imshow(np.rot90(dens_gaia), cmap=plt.cm.viridis, extent=[x.min(), x.max(), y_gaia.min(), y_gaia.max()], aspect='auto')
#plt.contour(X_gaia, Y_gaia, dens_gaia, colors='k', linewidth=0.01, levels = [10,30,50])
plt.text(-0.1, 10.0, 'Gaia DR1')
plt.xlim(-0.3,2.0)
#plt.ylim(-4.0,12.5)
plt.xlabel('(B-V)')
plt.yticks([])
plt.gca().invert_yaxis()

plt.subplots_adjust(wspace=0)
plt.savefig('Figure 3a-b kernel')
