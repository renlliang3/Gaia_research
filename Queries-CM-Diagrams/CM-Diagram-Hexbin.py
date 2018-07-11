import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.io.votable import parse
from astropy.table import Table
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

readresults = Table.read('Data/Updated-A-Query.fits',format='fits')
#readresults = Table.read('Reduced_DR2_Radial_Velocity_<10pc.csv',format='csv')
results = np.array(readresults)

x=results['bp_rp']#-results['e_bp_min_rp_val']
y_gaia=results['phot_g_mean_mag']+5*np.log10(results['parallax'])-10#-results['a_g_val']
k=np.linspace(-4, 1, 1000)

#counts_gaia,xbins_gaia,ybins_gaia,image_gaia = plt.hist2d(x,y_gaia,bins=250,normed=True,norm=LogNorm(),cmap = plt.cm.jet)
hb=plt.hexbin(x, y_gaia, extent=(0,0.33,-4,4), gridsize=150, bins='log', cmap = 'Blues')
plt.colorbar()
plt.vlines(x=0, ymin=-5, ymax=2,color=(1.0,0.2,0.3))
plt.vlines(x=0.33, ymin=-5, ymax=4,color=(1.0,0.2,0.3))
#plt.contour(counts_gaia.transpose(), extent=[xbins_gaia.min(),xbins_gaia.max(),ybins_gaia.min(),ybins_gaia.max()], colors='k', linewidth=0.01, levels = [0.1])
#plt.text(-0.1, 10.0, 'Gaia DR1')
plt.text(0.175, -3.5, 'A')
plt.plot(k,6*k+2,linestyle='-',color=(1.0,0.2,0.3))
plt.xlim(0,0.33)
plt.ylim(-4,4)
plt.xlabel(r'$(G_{BP}-G_{RP})$')
plt.ylabel(r'$M_G$')
plt.gca().invert_yaxis()

plt.savefig('CM-HexDiagram_Updated-A-Query')
