import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
import numpy as np

job = Gaia.launch_job_async("SELECT * \
FROM gaiadr2.gaia_source AS g \
WHERE g.parallax_over_error >= 5 \
AND g.parallax >= 2 \
AND g.bp_rp <= 0.9 \
AND g.phot_g_mean_mag + 5 * log10(g.parallax) - 10 <= 5", dump_to_file=True)
