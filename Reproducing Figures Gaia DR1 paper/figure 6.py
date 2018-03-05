import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
import numpy as np

job1 = Gaia.launch_job_async("select gaia.source_id, \
gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 as g_mag_abs, \
hip.b_v as b_min_v, \
sqrt(power(gaia.pmra,2)+power(gaia.pmdec,2))/gaia.parallax*4.74047 as vperp \
from gaiadr1.tgas_source as gaia \
inner join public.hipparcos_newreduction as hip \
on gaia.hip = hip.hip \
where gaia.parallax/gaia.parallax_error >= 5 and \
hip.e_b_v > 0.0 and hip.e_b_v <= 0.05 and \
2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux <= 0.05 and \
(gaia.parallax >= 10.0 or \
sqrt(power(gaia.pmra,2)+power(gaia.pmdec,2)) >= 200 or \
gaia.phot_g_mean_mag <= 7.5)", dump_to_file=True)

r1 = job1.get_results()

job2 = Gaia.launch_job_async("select gaia.source_id, \
gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 as g_mag_abs, \
0.85*(tycho2.bt_mag-tycho2.vt_mag) as b_min_v, \
sqrt(power(gaia.pmra,2)+power(gaia.pmdec,2))/gaia.parallax*4.74047 as vperp \
from gaiadr1.tgas_source as gaia \
inner join public.tycho2 as tycho2 \
on gaia.tycho2_id = tycho2.id \
where gaia.parallax/gaia.parallax_error >= 5 and \
sqrt(power(tycho2.e_bt_mag,2) + power(tycho2.e_vt_mag,2)) <= 0.05 and \
2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux <= 0.05 and \
(gaia.parallax >= 10.0 or \
sqrt(power(gaia.pmra,2)+power(gaia.pmdec,2)) >= 200 or \
gaia.phot_g_mean_mag <= 7.5)", dump_to_file=True)

r2 = job2.get_results()

job3 = Gaia.launch_job_async("select gaia.source_id, \
gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 as g_mag_abs, \
(urat.b_mag-urat.v_mag) as b_min_v, \
sqrt(power(gaia.pmra,2)+power(gaia.pmdec,2))/gaia.parallax*4.74047 as vperp \
from gaiadr1.tgas_source as gaia \
inner join public.tycho2 as tycho2 \
on gaia.tycho2_id = tycho2.id \
inner join gaiadr1.urat1_best_neighbour as uratxmatch \
on gaia.source_id = uratxmatch.source_id \
inner join gaiadr1.urat1_original_valid as urat \
on uratxmatch.urat1_oid = urat.urat1_oid \
where gaia.parallax/gaia.parallax_error >= 5 and \
sqrt(power(tycho2.e_bt_mag,2) + power(tycho2.e_vt_mag,2)) > 0.05 and \
sqrt(power(urat.b_mag_error,2) + power(urat.v_mag_error,2)) <= 0.05 and \
2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux <= 0.05 and \
(gaia.parallax >= 10.0 or \
sqrt(power(gaia.pmra,2)+power(gaia.pmdec,2)) >= 200 or \
gaia.phot_g_mean_mag <= 7.5)", dump_to_file=True)	

r3 = job3.get_results()

x1=r1['vperp']
color1=np.where(x1<10,'xkcd:dark purple',np.where(x1<20,'xkcd:sky blue',np.where(x1<50,'xkcd:turquoise',np.where(x1<100,'xkcd:deep green', np.where(x1<150,'xkcd:mustard', np.where(x1<200,'xkcd:orange','xkcd:magenta'))))))

#color1=np.where(x1>200,'xkcd:magenta',np.where(np.logical_and(150<=x1,x1<200),'xkcd:orange',np.where(np.logical_and(100<=x1,x1<150),'xkcd:mustard',np.where(np.logical_and(50<=x1,x1<100),'xkcd:deep green',np.where(np.logical_and(20<=x1,x1<50),'xkcd:turquoise',np.where(np.logical_and(10<=x1,x1<20),'xkcd:sky blue','xkcd:dark purple'))))))

x2=r2['vperp']
color2=np.where(x2<10,'xkcd:dark purple',np.where(x2<20,'xkcd:sky blue',np.where(x2<50,'xkcd:turquoise',np.where(x2<100,'xkcd:deep green', np.where(x2<150,'xkcd:mustard', np.where(x2<200,'xkcd:orange','xkcd:magenta'))))))

x3=r3['vperp']
color3=np.where(x3<10,'xkcd:dark purple',np.where(x3<20,'xkcd:sky blue',np.where(x3<50,'xkcd:turquoise',np.where(x3<100,'xkcd:deep green', np.where(x3<150,'xkcd:mustard', np.where(x3<200,'xkcd:orange','xkcd:magenta'))))))

plt.scatter(r1['b_min_v'], r1['g_mag_abs'], color=color1, s=0.1)
#plt.scatter(r2['b_min_v'], r2['g_mag_abs'], color=color2, s=0.1)
#plt.scatter(r3['b_min_v'], r3['g_mag_abs'], color=color3, s=0.1)
#plt.text(-0.1, 10.0, '')
plt.xlim(-0.3,2.0)
plt.ylim(-4.0,12.5)
plt.xlabel(r'$(B-V)$')
plt.ylabel(r'$M_G$')
plt.gca().invert_yaxis()

plt.savefig('Figure 6 hipparcos high vel first.eps', format='eps', dpi=1000)
