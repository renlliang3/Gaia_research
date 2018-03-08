import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
import numpy as np

r=[[],[],[]]

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

r[0] = job1.get_results()

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

r[1] = job2.get_results()

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

r[2] = job3.get_results()

def boundary(low,high,a):
	if type(low)!=str and type(high)!=str:
		return np.logical_and(low<=a,a<high)
	elif low=='None':
		return np.logical_and(a<high,0==0)
	elif high=='None':
		return np.logical_and(a>low,0==0)

x=[[],[],[]]
b_min_v=[[],[],[]]
g_mag_abs=[[],[],[]]

colors=['xkcd:sky blue','xkcd:teal','xkcd:deep green','xkcd:dark purple','xkcd:mustard','xkcd:orange','xkcd:magenta']
ranges=[[10,20],[20,50],[50,100],['None',10],[100,150],[150,200],[200,'None']]

for i in range(len(colors)):
	for j in range(len(x)):
		x[j]=r[j]['vperp']
		b_min_v[j].append(r[j]['b_min_v'][np.where(boundary(ranges[i][0],ranges[i][1],x[j]))])
		g_mag_abs[j].append(r[j]['g_mag_abs'][np.where(boundary(ranges[i][0],ranges[i][1],x[j]))])

for i in range(len(colors)):
	for j in range(len(x)):
		plt.scatter(b_min_v[j][i],g_mag_abs[j][i],color=colors[i],s=0.00001)

plt.xlim(-0.3,2.0)
plt.ylim(-4.0,12.5)
plt.xlabel(r'$(B-V)$')
plt.ylabel(r'$M_G$')
plt.gca().invert_yaxis()

plt.savefig('Figure 6 plotting order from low to high velocities (lowest fourth)', format='eps', dpi=1000)
