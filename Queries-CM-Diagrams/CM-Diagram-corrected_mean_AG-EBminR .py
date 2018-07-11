import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.io.votable import parse
from astropy.table import Table
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

extra_data = np.genfromtxt("Data/Mean_extinction_excess_SourceID_500pc_AG-4.4_BminR-1.7.dat", names = True, dtype=None)
source_id = extra_data['source_id']

readresults = Table.read('Data/All-star-500pc-AG-4.4-BminR-1.7.fits',format='fits')
results = np.array(readresults)

#Find the nonmatching source_id
nonmatches = np.array([])
j=0

for i in range(len(source_id)):
	not_found = True
	while not_found:
		if source_id[i]!=results['source_id'][j]:
			nonmatches = np.append(nonmatches,j)
			j+=1
		else:
			not_found = False
			j+=1
	print(j)

#Delete the rows with source_id that have no AG and EBminR
for k in range(len(nonmatches)):
	results = np.delete(results,nonmatches[k]-k)

#Check if the deletion was succesful
j=0
nonmatches_check = np.array([])
for i in range(len(source_id)):
	not_found = True
	while not_found:
		if source_id[i]!=results['source_id'][j]:
			nonmatches_check = np.append(nonmatches_check,j)
			j+=1
		else:
			not_found = False
			j+=1
	print(j)

x=results['bp_rp']-extra_data['EBminR']
y_gaia=results['phot_g_mean_mag']+5*np.log10(results['parallax'])-10-extra_data['AG']
k=np.linspace(-4, 1, 1000)

counts_gaia,xbins_gaia,ybins_gaia,image_gaia = plt.hist2d(x,y_gaia,bins=600,normed=True,norm=LogNorm(),cmap = plt.cm.jet)
plt.colorbar()
plt.vlines(x=0, ymin=-5, ymax=2,color=(1.0,0.2,0.3))
plt.vlines(x=0.33, ymin=-5, ymax=4,color=(1.0,0.2,0.3))
#plt.contour(counts_gaia.transpose(), extent=[xbins_gaia.min(),xbins_gaia.max(),ybins_gaia.min(),ybins_gaia.max()], colors='k', linewidth=0.01, levels = [0.1])
plt.text(-0.6, -3.5, 'OB')
plt.text(0.1, -3.5, 'A')
plt.plot(k,3*k+2.1,linestyle='-',color=(1.0,0.2,0.3))
plt.xlim(np.min(x),1.7)
plt.ylim(-4,4.4)
plt.xlabel(r'$(G_{BP}-G_{RP})-<E(B-R)>$')
plt.ylabel(r'$M_G-<A_G>$')
plt.gca().invert_yaxis()

plt.savefig('CM-Diagram_All-star-500pc-AG-4.4-BminR-1.7-Corrected_mean_AG_EBminR.png')
