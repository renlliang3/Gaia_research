import numpy as np
import astropy.coordinates as asc
from astropy import units as u

gammavelorum_icrs = asc.SkyCoord(ra='08h09m31.95013s', dec=(-47-(20/60)-(11.7108/3600))*u.degree, distance=336*u.pc, frame='icrs')
gammavelorum_galactic = gammavelorum_icrs.galactic
gammavelorum_cartesian_galactic = asc.spherical_to_cartesian(gammavelorum_galactic.distance,gammavelorum_galactic.b,gammavelorum_galactic.l)
print("gamma Velorum")
print(gammavelorum_cartesian_galactic)
print("")

#velaOB2_galactic = asc.SkyCoord(l=263*u.degree, b=-7*u.degree, distance=450*u.pc, frame='galactic')
velaOB2_galactic = asc.SkyCoord(l=263*u.degree, b=-7*u.degree, distance=410*u.pc, frame='galactic')
velaOB2_cartesian_galactic = asc.spherical_to_cartesian(velaOB2_galactic.distance,velaOB2_galactic.b,velaOB2_galactic.l)
print("Vela OB2")
print(velaOB2_cartesian_galactic)
print("")

#NGC2547_ICRS = asc.SkyCoord(ra=360*(8/24 + 10.7/(24*60))*u.degree,dec=-(49+16/60)*u.degree, distance=((1.97*10**3)/3.26156)*u.pc, frame='icrs')
NGC2547_ICRS = asc.SkyCoord(ra=360*(8/24 + 10.7/(24*60))*u.degree,dec=-(49+16/60)*u.degree, distance=410*u.pc, frame='icrs')
NGC2547_galactic = NGC2547_ICRS.galactic
NGC2547_cartesian_galactic = asc.spherical_to_cartesian(NGC2547_galactic.distance,NGC2547_galactic.b,NGC2547_galactic.l)
print("NGC2547")
print(NGC2547_cartesian_galactic)
print("")

#trumpler10_ICRS = asc.SkyCoord(ra=360*(8/24+47.8/(24*60))*u.degree, dec=-(42+29/60)*u.degree, distance=337.25*u.pc, frame='icrs')
trumpler10_ICRS = asc.SkyCoord(ra=360*(8/24+47.8/(24*60))*u.degree, dec=-(42+29/60)*u.degree, distance=366*u.pc, frame='icrs')
trumpler10_galactic = trumpler10_ICRS.galactic
trumpler10_cartesian_galactic = asc.spherical_to_cartesian(trumpler10_galactic.distance,trumpler10_galactic.b,trumpler10_galactic.l)
print("Trumpler 10")
print(trumpler10_cartesian_galactic)
print("")

#NGC2516_ICRS = asc.SkyCoord(ra=360*(7/24+58/(24*60)+20/(24*3600))*u.degree, dec=-(60+52/60)*u.degree, distance=(1300/3.26156)*u.pc, frame='icrs')
NGC2516_ICRS = asc.SkyCoord(ra=360*(7/24+58/(24*60)+20/(24*3600))*u.degree, dec=-(60+52/60)*u.degree, distance=373*u.pc, frame='icrs')
NGC2516_galactic = NGC2516_ICRS.galactic
NGC2516_cartesian_galactic = asc.spherical_to_cartesian(NGC2516_galactic.distance,NGC2516_galactic.b,NGC2516_galactic.l)
print("NGC2516")
print(NGC2516_cartesian_galactic)
print("")

test = asc.SkyCoord(ra=360*(7/24+58/(24*60)+20/(24*3600))*u.degree, dec=-(60+52/60)*u.degree, distance=373*u.pc, frame='icrs', obstime='J2015')
test_galactic = test.galactic
test_cartesian_galactic = asc.spherical_to_cartesian(test_galactic.distance,test_galactic.b,test_galactic.l)
print("test obstime J2015")
print(test_cartesian_galactic)
print("")

x=np.array([[360*(7/24+58/(24*60)+20/(24*3600)), 360*(8/24+47.8/(24*60))], [-(60+52/60), -(42+29/60)], [373, 366]])
test2 = asc.SkyCoord(ra=x[0]*u.degree, dec=x[1]*u.degree, distance=x[2]*u.pc, frame='icrs', obstime='J2015')
test2_galactic = test2.galactic
test2_cartesian_galactic = asc.spherical_to_cartesian(test2_galactic.distance,test2_galactic.b,test2_galactic.l)
print("test array input (NGC2516,trumpler10)")
print(test2_cartesian_galactic)
