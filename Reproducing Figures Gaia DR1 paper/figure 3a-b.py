import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
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
print(r['g_mag_abs_gaia'])



plt.subplot(1, 2, 1)
plt.scatter(r['b_v'], r['g_mag_abs_hip'], color='g', s=0.1)
plt.text(-0.1, 10.0, 'Hipparcos')
plt.xlim(-0.3,2.0)
plt.ylim(-4.0,12.5)
plt.xlabel('(B-V)')
plt.ylabel(r'$M_G$')
plt.gca().invert_yaxis()

plt.subplot(1, 2, 2)
plt.scatter(r['b_v'], r['g_mag_abs_gaia'], color='g', s=0.1)
plt.text(-0.1, 10.0, 'Gaia DR1')
plt.xlim(-0.3,2.0)
plt.ylim(-4.0,12.5)
plt.xlabel('(B-V)')
plt.yticks([])
plt.gca().invert_yaxis()

plt.subplots_adjust(wspace=0)
plt.savefig('Figure 3a-b')
