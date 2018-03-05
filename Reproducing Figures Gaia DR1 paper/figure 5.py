import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
import numpy as np

job = Gaia.launch_job_async("select gaia.source_id, \
gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 as g_mag_abs, \
gaia.phot_g_mean_mag-tmass.ks_m as g_min_ks \
from gaiadr1.tgas_source as gaia \
inner join gaiadr1.tmass_best_neighbour as xmatch \
on gaia.source_id = xmatch.source_id \
inner join gaiadr1.tmass_original_valid as tmass \
on tmass.tmass_oid = xmatch.tmass_oid \
where gaia.parallax/gaia.parallax_error >= 5 and ph_qual = ’AAA’ and \
sqrt(power(2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux,2)) <= 0.05 and \
sqrt(power(2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux,2) \
+ power(tmass.ks_msigcom,2)) <= 0.05", dump_to_file=True)

print(job)

r = job.get_results()

plt.scatter(r['g_min_ks'], r['g_mag_abs'], color='b', s=0.01)
plt.text(-0.1, 10.0, 'Gaia DR1')
plt.xlim(-0.6,4.0)
plt.ylim(-3.0,12.5)
plt.xlabel(r'$(G-K_S)$')
plt.ylabel(r'$M_G$')
plt.gca().invert_yaxis()

plt.savefig('Figure 5 py')
