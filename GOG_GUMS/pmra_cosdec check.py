import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.table import Table
from scipy import stats
from mayavi import mlab
import astropy.coordinates as asc
import numpy as np
import random

readresults = Table.read('GOG_GUMS_DR2_7e6.fits',format='fits')
results = np.array(readresults)

dist = results['dist']
list = []

for i in range(len(dist)):
	if (dist[i] < 650):
		list.append(i)

list = np.array(list)

dist = np.array([dist[i] for i in list])

source_id = np.array([results['source_id'][i] for i in list])
ra = np.array([results['ra'][i] for i in list])
dec = np.array([results['dec'][i] for i in list])
varpi = np.array([results['varpi'][i] for i in list])
e_varpi = np.array([results['e_varpi'][i] for i in list])
pmra = np.array([results['pmra'][i] for i in list])
e_pmra = np.array([results['e_pmra'][i] for i in list])
pmdec = np.array([results['pmdec'][i] for i in list])
e_pmdec = np.array([results['e_pmdec'][i] for i in list])
vrad = ((60/22)**0.5)*np.array([results['vrad'][i] for i in list])
e_vrad = np.array([results['e_vrad'][i] for i in list])
G_mag = np.array([results['G_Mag'][i] for i in list])
GRVS_Mag = np.array([results['GRVS_Mag'][i] for i in list])
GBP_Mag = np.array([results['GBP_Mag'][i] for i in list])
GRP_Mag = np.array([results['GRP_Mag'][i] for i in list])
dist = np.array([results['dist'][i] for i in list])

pmra_cosdec = pmra*np.cos(dec*np.pi/180)

indices=[]
for i in range(len(dec)):
        if (np.abs(dec[i])>80):
                indices.append(i)

highdecpmra=np.array([pmra[i] for i in indices])
highdecpmra_cosdec=np.array([pmra_cosdec[i] for i in indices])

indices=[]
for i in range(len(dec)):
        if (40<np.abs(dec[i])<50):
                indices.append(i)
                
mediumdecpmra=np.array([pmra[i] for i in indices])
mediumdecpmra_cosdec=np.array([pmra_cosdec[i] for i in indices])

indices=[]
for i in range(len(dec)):
        if (np.abs(dec[i])<10):
                indices.append(i)
lowdecpmra=np.array([pmra[i] for i in indices])
lowdecpmra_cosdec=np.array([pmra_cosdec[i] for i in indices])

print('Variance of pmra for 80<|dec|<90: '+str(np.var(highdecpmra)))
print('Variance of pmra_cosdec for 80<|dec|<90: '+str(np.var(highdecpmra_cosdec)))
print('Variance of pmra for 40<|dec|<50: '+str(np.var(mediumdecpmra)))
print('Variance of pmra_cosdec for 40<|dec|<50: '+str(np.var(mediumdecpmra_cosdec)))
print('Variance of pmra for |dec|<10: '+str(np.var(lowdecpmra)))
print('Variance of pmra_cosdec for |dec|<10: '+str(np.var(lowdecpmra_cosdec)))
