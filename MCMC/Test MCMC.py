import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op
import corner
from scipy.stats import reciprocal

# Generate synthetic survey for N = 400 stars with absolute magnitudes drawn from Gaussian luminosity distribution function with mu_M=9 and variance_M=0.49
# and parallaxes drawn from a uniform distrubituion between d=1 and d=100pc
N = 400
mu_Mag_true=9
var_Mag_true=0.49
sigma_Mag_true=var_Mag_true**0.5
np.random.seed(42)
Magnitudes_abs_true = np.random.normal(mu_Mag_true,np.sqrt(var_Mag_true),N)
parallaxes_true = reciprocal.rvs(0.01,1,size=N)
magnitudes_app_true = Magnitudes_abs_true - 5*np.log10(parallaxes_true) - 5

#Simulate synthetic observations for each of the stars, we do this by drawing from a Gaussian around the true mean we have generated in our synthetic survey for the True values of the stars, with variance defined in the functions below.
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

sigma_parallaxes_obs = [sigma_varpi(x) for x in magnitudes_app_true]
sigma_magnitudes_app_obs = [sigma_mag(x) for x in magnitudes_app_true]

magnitudes_app_obs = [np.random.normal(x,sigma_mag(x)) for x in magnitudes_app_true]
parallaxes_obs=[]
for i in range(len(parallaxes_true)):
	parallaxes_obs.append(np.random.normal(parallaxes_true[i],sigma_parallaxes_obs[i]))

#The priors of our Bayes setup model. Which is a Gaussian with something added 
def lnprob(x, mu, sigma):
	y = -0.5*(((x-mu)/sigma)**2)
	if np.isnan(y):
		return -np.inf
	else:
		return y

def lnpriormu_M(mu_M):
	if -5 < mu_M < 16:
		return np.log(1)
	else:
		return -np.inf #np.log(0)

def lnpriorsigma_M(sigma_M):
	if 0.1 < sigma_M < 1.0:
		return np.log(1/(sigma_M**2))
	else:
		return -np.inf #np.log(0)

#For these two priors we will need them for each star so 400 different priors multiplied together for each of prior types
def lnpriorparallax_true(parallax):
	if 10**-3 < parallax < 1:
		return np.log(parallax**-4)
	else:
		return -np.inf

def lnpriorMagnitude_abs_true(M,mu_M,sigma_M):
	return -np.log(sigma_M)+lnprob(M,mu_M,sigma_M)

#The likelihoods of our bayes setup, again 400 different ones multiplied together
def lnlikelihoodparallax_obs(parallax_obs,parallax_true,sigma_parallax_obs):
	return lnprob(parallax_obs,parallax_true,sigma_parallax_obs)

#print(lnlikelihoodparallax_obs(5,10,3))

def lnlikelihoodmagnitude_app_obs(magnitude_app_obs,magnitude_app_true,sigma_magnitude_app_obs):
	return lnprob(magnitude_app_obs,magnitude_app_true,sigma_magnitude_app_obs)

#The posterior of our Bayes setup. Order of input variables needs to be parameters, observations, which it is.
def lnposterior(theta,magnitudes_app_obs,parallaxes_obs,sigma_magnitudes_app_obs,sigma_parallaxes_obs):
	muMTrue, varMTrue, pTrue, MTrue = theta[0], theta[1], theta[2:N+2], theta[N+2:]
	#print(muMTrue)	
	#print(varMTrue)
	mappTrue = MTrue - 5*np.log10(pTrue) - 5
	#print(np.shape(sigma_parallaxes_obs))

	print(mappTrue)

	#for i in range(len(parallaxes_obs)):
		#print(np.log(parallaxes_obs[i]**-4))

	#for i in range(len(mappTrue)):
	#	if np.isnan(mappTrue[i]):
	#		mappTrue[i]=-9.7

	holder=0

	holder1=0
	holder2=0
	holder3=0
	holder4=0
	
	print(np.min(pTrue))
	print(np.max(pTrue))
	
	holder+=-4.0*np.sum(np.log(pTrue*1000))
	#print(magnitudes_app_obs[397])
	#print(mappTrue[397])
	#print(sigma_magnitudes_app_obs[397])
	for i in range(N):
		holder+=lnlikelihoodmagnitude_app_obs(magnitudes_app_obs[i],mappTrue[i],sigma_magnitudes_app_obs[i])
		#print(lnlikelihoodmagnitude_app_obs(magnitudes_app_obs[i],mappTrue[i],sigma_magnitudes_app_obs[i]))
		holder+=lnlikelihoodparallax_obs(parallaxes_obs[i],pTrue[i],sigma_parallaxes_obs[i])
		#print(lnlikelihoodparallax_obs(parallaxes_obs[i],pTrue[i],sigma_parallaxes_obs[i]))
		holder+=lnpriorparallax_true(pTrue[i])
		#holder3+=lnpriorparallax_true(pTrue[i])
		#print(lnpriorparallax_true(pTrue[i]))
		holder+=lnpriorMagnitude_abs_true(MTrue[i],muMTrue,varMTrue**0.5)
	#print(holder1)
	#print(holder2)
	print(holder3)
	#print(holder4)
	holder+=lnpriormu_M(muMTrue)
	#print(muMTrue)
	#print(lnpriormu_M(muMTrue))
	holder+=lnpriorsigma_M(varMTrue**0.5)
	#print(varMTrue**0.5)
	#print(lnpriorsigma_M(varMTrue**0.5))
	return holder

def UniformSpaceDensityGaussianLFBookemcee(posWalker, posteriorDict, observations, observationalErrors):
	lnPosterior=0.0
	meanAbsMag=posWalker[0]
	if (meanAbsMag<-5 or meanAbsMag>16):
		return -np.inf
	#lnPosterior=lnPosterior-np.log(posteriorDict['muHigh']-posteriorDict['muLow'])

	variance=posWalker[1]
	if (variance<0.01 or variance>1):
		return -np.inf
	lnPosterior=lnPosterior - np.log(variance)

	parallaxes=posWalker[2:numStarsInSurvey+2]
	if (np.any(parallaxes < posteriorDict['minParallax']) or np.any(parallaxes > posteriorDict['maxParallax'])):
		return -np.inf
	lnPosterior = lnPosterior - 4.0*np.sum(np.log(parallaxes))

	absoluteMagnitudes=posWalker[numStarsInSurvey+2:]
	apparentMagnitudes=absoluteMagnitudes-5.0*np.log10(parallaxes)+10.0

	y=np.concatenate((parallaxes,apparentMagnitudes,absoluteMagnitudes))
	inv_var=np.diagflat(np.concatenate((observationalErrors, 1.0/np.repeat(variance,numStarsInSurvey))))
	means=np.concatenate((observations, np.repeat(meanAbsMag,numStarsInSurvey)))
	diff=y-means
	lnPosterior = lnPosterior - 0.5*np.dot(diff,np.dot(inv_var,diff))
	lnPosterior = lnPosterior - 0.5*numStarsInSurvey*np.log(variance)
	return lnPosterior

#lnposterior([10 for i in range(802)],[5 for i in range(802)],

#Setting up the MCMC method, starting position will be an educated guess from the observations we have.
ndim = 2*N+2
nwalkers = ndim * 10

#Initial guess for the parallaxes, we clip the parallaxes between these values because they are supposed to be between these values. This will solve any NaN problems aswell
initialpTrue = np.clip(parallaxes_obs,10**-3, 1)

#Some negative parallaxes cause NaN. Replace the values with -15. (-15 was smallest number I found in the list)
initialMTrue = magnitudes_app_obs + 5*np.log10(initialpTrue) + 5
#for i in range(len(initialMTrue)):
#	if np.isnan(initialMTrue[i]):
#		initialMTrue[i]=-15

initialmuMTrue = np.mean(initialMTrue)
initialvarMTrue = (1-0.01)/(np.log(1)-np.log(0.01)) #np.var(initialMTrue)

#parameter order is: [mu_MTrue,var_MTrue,parallaxes_true[0],...,parallaxes_true[399],Mag_true[0],...,Mag_true[399]]
initialParameters = np.concatenate((np.array([initialmuMTrue, initialvarMTrue]),initialpTrue, initialMTrue))

#Initialize the walkers around this initialguess, but each in a slightly different position. This does deviate from the initialization position that Anthony uses.
initialPos = [initialParameters + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

args=(magnitudes_app_obs,parallaxes_obs,sigma_magnitudes_app_obs,sigma_parallaxes_obs)
for i in range(1):
	print(lnposterior(initialPos[i],*args))
