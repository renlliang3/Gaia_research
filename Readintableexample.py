from astropy.table import Table

r=Table.read('Fig5-result.vot',format='votable')

print(len(r['g_min_ks']))
