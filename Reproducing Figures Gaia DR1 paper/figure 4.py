import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import numpy as np

job = Gaia.launch_job_async("select gaia.source_id, gaia.hip, \
gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 as g_mag_abs_gaia, \
gaia.phot_g_mean_mag+5*log10(hip.plx)-10 as g_mag_abs_hip \
from gaiadr1.tgas_source as gaia \
inner join public.hipparcos_newreduction as hip \
on gaia.hip = hip.hip \
where gaia.parallax/gaia.parallax_error >= 5 and \
hip.plx/hip.e_plx >= 5 and \
hip.e_b_v > 0.0 and hip.e_b_v <= 0.05 and \
hip.b_v >= 1.0 and hip.b_v <= 1.1 and \
2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux <= 0.05", dump_to_file=True)

print(job)

r = job.get_results()

mag_gaia=r['g_mag_abs_gaia']
mag_hip=r['g_mag_abs_hip']
Xaxis_plot=np.linspace(-3,8.5,1000)

kde = KernelDensity(bandwidth=0.2, kernel='epanechnikov')
log_dens_mag_gaia = kde.fit(mag_gaia[:,None]).score_samples(Xaxis_plot[:,None])
log_dens_mag_hip = kde.fit(mag_hip[:,None]).score_samples(Xaxis_plot[:,None])

plt.plot(Xaxis_plot, np.exp(log_dens_mag_gaia), 'g', label='Gaia DR1')
plt.plot(Xaxis_plot, np.exp(log_dens_mag_hip), 'b', label='Hipparcos')
plt.legend(loc='upper right')
plt.ylim(0,1.1)
plt.xlabel(r'$M_G$')
plt.ylabel('Density')
plt.savefig('Figure 4')
