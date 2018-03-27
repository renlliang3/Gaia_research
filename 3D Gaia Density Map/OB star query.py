import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

job = Gaia.launch_job_async("select gaia.*, \
gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 as g_mag_abs \
from gaiadr1.gaia_source as gaia \
where gaia.parallax_error/gaia.parallax <= 0.2 and \
gaia.parallax >= 1.538 and gaia.phot_g_mean_mag+5*log10(gaia.parallax)-10 < 1.7", dump_to_file = True)
#2.5/log(10)*gaia.phot_g_mean_flux_error/gaia.phot_g_mean_flux <= 0.05", dump_to_file=True
