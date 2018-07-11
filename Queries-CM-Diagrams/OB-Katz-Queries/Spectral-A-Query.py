import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

job = Gaia.launch_job_async("SELECT * \
FROM gaiadr2.gaia_source AS g, gaiadr2.tmass_best_neighbour AS tbest, gaiadr1.tmass_original_valid AS tmass \
WHERE g.source_id = tbest.source_id AND tbest.tmass_oid = tmass.tmass_oid \
AND tmass.ph_qual = 'AAA' AND tmass.j_msigcom <= 0.05 AND tmass.h_msigcom <= 0.05 AND tmass.ks_msigcom <= 0.05 \
AND g.parallax >= 2 \
AND g.parallax_over_error >= 5 \
AND g.bp_rp - g.e_bp_min_rp_val >= -0.04 \
AND g.bp_rp - g.e_bp_min_rp_val <= 0.33 \
AND g.phot_g_mean_mag + 5 * log10(g.parallax) - 10 - g.a_g_val >= 1.1 \
AND g.phot_g_mean_mag + 5 * log10(g.parallax) - 10 - g.a_g_val <= 2.3 \
AND tmass.j_m - tmass.h_m <= 0.14 * (g.phot_g_mean_mag - tmass.ks_m) + 0.02 \
AND tmass.j_m - tmass.ks_m <= 0.23 * (g.phot_g_mean_mag - tmass.ks_m)", dump_to_file=True)
