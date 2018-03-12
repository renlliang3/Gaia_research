import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op
import corner
from scipy.stats import reciprocal

# Generate synthetic survey for N = 400 stars with absolute magnitudes drawn from Gaussian luminosity distribution function with mu_M=9 and variance_M=0.49
# and parallaxes drawn from a uniform distrubituion between d=1 and d=100pc
def sigma_varpi(m_i):
	m_0=5
	a_varpi=0.2*10**-3
	b_varpi=0.2*10**-3
	if m_i >= m_0:
		return a_varpi*10**(0.2*(m_i-m_0))
	else:
		return b_varpi

def sigma_mag(m_i):
	m_0=5
	a_m=0.006
	b_m=0.001
	if m_i >= m_0:
		return a_m*10**(0.2*(m_i-m_0))
	else:
		return b_m

N = 400
mu_M=9
var_M=0.49
Magnitudes = np.random.normal(mu_M,np.sqrt(var_M),400)
Parallaxes = reciprocal.rvs(0.01,1,size=N)
sigma_parallaxes = [sigma_varpi(x) for x in Magnitudes] 
sigma_magnitudes = [sigma_mag(x) for x in Magnitudes]

#The likelihood of our Bayes setup model. Which is a Gaussian with something added 

def lnprob(x, mu, sigma):
    diff = x-mu
    return -(((x-mu)/sigma)**2)/2.0

def lnpriormu_M(mu_M):
	if -5 < mu_m < 16:
		return 1/21.
	else:
		return 0

def lnpriorsigma_M(sigma_M):
	if 0.1 < sigma_M < 1.0:
		return (1/9)*(1/sigma_m**2)
	else:
		return 0


#Still need to simulate the observations...
likelihood=1
for i in range(N):
	likelihood*=np.exp(-(((varpi_obs-parallaxes[i])/sigma_parallaxes[i])**2)/2.0)*np.exp(-(((m_i-Magnitudes[i]+5.*np.log10(varpi_obs)+5))/sigma_magnitudes[i])**2)/2.0)



plt.subplot(1,4,1)
plt.hist(Magnitudes,density=True)
plt.subplot(1,4,2)
plt.hist(Parallaxes,density=True)
plt.subplot(1,4,3)
plt.hist(sigma_parallaxes,density=True)
plt.subplot(1,4,4)
plt.hist(sigma_magnitudes,density=True)
plt.show()
