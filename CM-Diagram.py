import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.io.votable import parse
from astropy.table import Table
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

readresults = Table.read('Data/OBvrad1000-Katz.vot',format='votable')
#readresults = Table.read('Reduced_DR2_Radial_Velocity_<10pc.csv',format='csv')
results = np.array(readresults)

x=results['bp_rp']#-results['e_bp_min_rp_val']
y_gaia=results['phot_g_mean_mag']+5*np.log10(results['parallax'])-10#-results['a_g_val']

counts_gaia,xbins_gaia,ybins_gaia,image_gaia = plt.hist2d(x,y_gaia,bins=250,normed=True
                                      ,norm=LogNorm()
                                      , cmap = plt.cm.viridis)
plt.colorbar()
#plt.contour(counts_gaia.transpose(), extent=[xbins_gaia.min(),xbins_gaia.max(),ybins_gaia.min(),ybins_gaia.max()], colors='k', linewidth=0.01, levels = [0.1])
#plt.text(-0.1, 10.0, 'Gaia DR1')
plt.xlim(x.min(),x.max())
plt.ylim(y_gaia.min(),y_gaia.max())
plt.xlabel(r'$(G_{BP}-G_{RP})-E(G_{BP}-G_{RP})$')
plt.ylabel(r'$M_G-A_G$')
plt.gca().invert_yaxis()

plt.savefig('CM-Diagram_OBvrad1000-Katz-corrected_colour_extinction')
