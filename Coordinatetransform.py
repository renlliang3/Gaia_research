import numpy as np
import astropy.coordinates as asc
from astropy import units as u

gammavelorum_icrs = asc.SkyCoord(ra='08h09m31.95013s', dec=(-47-(20/60)-(11.7108/3600))*u.degree, distance=336*u.pc, frame='icrs')
print(gammavelorum_icrs)
gammavelorum_galactic = gammavelorum_icrs.galactic
print(gammavelorum_galactic)
gammavelorum_cartesian_galactic = asc.spherical_to_cartesian(gammavelorum_galactic.distance,gammavelorum_galactic.b,gammavelorum_galactic.l)
print(gammavelorum_cartesian_galactic)

#velaOB2_galactic = asc.SkyCoord(l=263*u.degree, b=-7*u.degree, distance=450*u.pc, frame='galactic')
velaOB2_galactic = asc.SkyCoord(l=263*u.degree, b=-7*u.degree, distance=410*u.pc, frame='galactic')
print(velaOB2_galactic)
velaOB2_cartesian_galactic = asc.spherical_to_cartesian(velaOB2_galactic.distance,velaOB2_galactic.b,velaOB2_galactic.l)
print(velaOB2_cartesian_galactic)

#NGC2547_ICRS = asc.SkyCoord(ra=360*(8/24 + 10.7/(24*60))*u.degree,dec=-(49+16/60)*u.degree, distance=((1.97*10**3)/3.26156)*u.pc, frame='icrs')
NGC2547_ICRS = asc.SkyCoord(ra=360*(8/24 + 10.7/(24*60))*u.degree,dec=-(49+16/60)*u.degree, distance=410*u.pc, frame='icrs')
NGC2547_galactic = NGC2547_ICRS.galactic
print(NGC2547_galactic)
NGC2547_cartesian_galactic = asc.spherical_to_cartesian(NGC2547_galactic.distance,NGC2547_galactic.b,NGC2547_galactic.l)
print(NGC2547_cartesian_galactic)
