import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.io.votable import parse
from astropy.table import Table
from astropy.io import fits
import numpy as np

"""
readresults = Table.read('Data/All-star-500pc-AG-4.4-BminR-1.7-crossmatched.fits',format='fits')
extra_data = Table.read("Data/Mean_extinction_excess_SourceID_500pc_AG-4.4_BminR-1.7.fits",format='fits')
"""

"""
t = Table.read("Data/Mean_extinction_excess_SourceID_500pc_AG-4.4_BminR-1.7.fits",format='fits')
t.remove_row(0)
t.rename_column('col1', 'source_id_2')
t.rename_column('col2', 'AG')
t.rename_column('col3', 'EBminR')
t.write('Data/Mean_extinction_excess_SourceID_500pc_AG-4.4_BminR-1.7.fits', format='fits')
"""

with fits.open('Data/All-star-500pc-AG-4.4-BminR-1.7-crossmatched.fits') as hdul1:
    with fits.open('Data/Mean_extinction_excess_SourceID_500pc_AG-4.4_BminR-1.7.fits') as hdul2:
        
        new_columns = hdul1[1].columns + hdul2[1].columns
        new_hdu = fits.BinTableHDU.from_columns(new_columns)


new_hdu.writeto('Concatenated.fits')
