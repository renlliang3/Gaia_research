import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.io.votable import parse
from astroquery.gaia import Gaia
from matplotlib.colors import LogNorm
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
where gaia.parallax/gaia.parallax_error >= 5 and ph_qual = 'AAA' and \
sqrt(power(2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux,2)) <= 0.05 and \
sqrt(power(2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux,2) \
+ power(tmass.ks_msigcom,2)) <= 0.05", dump_to_file=True)

r = job.get_results()

x=r['g_min_ks']
y_gaia=r['g_mag_abs']

counts_gaia,xbins_gaia,ybins_gaia,image_gaia = plt.hist2d(x,y_gaia,bins=250
                                      ,norm=LogNorm()
                                      , cmap = plt.cm.viridis)
plt.contour(counts_gaia.transpose(), extent=[xbins_gaia.min(),xbins_gaia.max(),ybins_gaia.min(),ybins_gaia.max()], colors='k', linewidth=0.01, levels = [10,30,50,70,90])
plt.text(-0.1, 10.0, 'Gaia DR1')
plt.xlim(-0.6,4.0)
plt.ylim(-3.0,12.5)
plt.xlabel(r'$(G-K_S)$')
plt.ylabel(r'$M_G$')
plt.gca().invert_yaxis()

plt.savefig('Figure 5')
