from astropy.table import Table
from numpy as np

r=Table.read('Fig5-result.vot',format='votable')

print(len(r['g_min_ks']))
