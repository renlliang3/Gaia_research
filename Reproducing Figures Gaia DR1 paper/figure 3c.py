import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

job = Gaia.launch_job_async("select gaia.source_id, gaia.hip, \
gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 as g_mag_abs, \
hip.b_v \
from gaiadr1.tgas_source as gaia \
inner join public.hipparcos_newreduction as hip \
on gaia.hip = hip.hip \
where gaia.parallax/gaia.parallax_error >= 5 and \
hip.e_b_v > 0.0 and hip.e_b_v <= 0.05 and \
2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux <= 0.05", dump_to_file=True)

r = job.get_results()
x=r['b_v']
y_gaia=r['g_mag_abs']

counts_gaia,xbins_gaia,ybins_gaia,image_gaia = plt.hist2d(x,y_gaia,bins=300,normed=True
                                      ,norm=LogNorm()
                                      , cmap = plt.cm.viridis)
plt.colorbar()
plt.contour(counts_gaia.transpose(), extent=[xbins_gaia.min(),xbins_gaia.max(),ybins_gaia.min(),ybins_gaia.max()], colors='k', linewidth=0.01, levels = [0.1,0.3,0.5])
plt.text(-0.1, 10.0, 'Gaia DR1')
plt.xlim(-0.3,2.0)
plt.ylim(-4.0,12.5)
plt.xlabel('(B-V)')
plt.gca().invert_yaxis()

plt.savefig('Figure 3c')
