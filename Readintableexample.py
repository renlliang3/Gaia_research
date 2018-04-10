from astropy.table import Table
import numpy as np

r = Table.read('Fig5-result.vot',format='votable')
r = np.array(r)
