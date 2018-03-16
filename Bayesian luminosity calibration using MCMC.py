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
    return -(((x-mu)/sigma)**2)/2.0

def lnpriormu_M(mu_M):
	if -5 < mu_m < 16:
		return np.log(1/21.)
	else:
		return -np.inf #np.log(0)

def lnpriorsigma_M(sigma_M):
	if 0.1 < sigma_M < 1.0:
		return np.log((1/9)*(1/sigma_m**2))
	else:
		return -np.inf #np.log(0)

#For these two priors we will need them for each star so 400 different priors multiplied together for each of prior types
def lnpriorparallax_true(parallax):
	return np.log(parallax**-4)

def lnpriorMagnitude_abs_true(M,mu_M,sigma_M):
	return -np.log(sigma_M*np.sqrt(2*np.pi))+lnprob(M,mu_m,sigma_M)

#The likelihoods of our bayes setup, again 400 different ones multiplied together
def lnlikelihoodparallax_obs(parralax_obs,parallax_true,sigma_parallax_obs):
	return lnprob(parallax_obs,parallax_true,sigma_parallax_obs)

print(lnlikelihoodparallax_obs(5,10,3))

def lnlikelihoodmagnitude_app_obs(magnitude_app_obs,magnitude_app_true,sigma_magnitude_app_obs):
	return lnprob(magnitude_app_obs,magnitude_app_true,sigma_magnitude_app_obs)

#The posterior of our Bayes setup. Order of input variables needs to be parameters, observations, which it is.
def lnposterior(theta,magnitudes_app_obs,parallaxes_obs,sigma_magnitudes_app_obs,sigma_parallaxes_obs):
	muMTrue, varMTrue, MTrue, pTrue = theta[0], theta[1], theta[2:402], theta[402:802]	
	mappTrue = MTrue - 5*np.log10(pTrue) - 5

	for i in range(len(mappTrue)):
		if np.isnan(mappTrue[i]):
			initialMTrue[i]=3

	holder=0

	print(len(magnitudes_app_obs))
	print(len(mappTrue))
	print(len(sigma_magnitudes_app_obs))
	for i in range(N):
		holder+=lnlikelihoodmagnitude_app_obs(magnitudes_app_obs[i],mappTrue[i],sigma_magnitudes_app_obs[i])
		holder+=lnlikelihoodparallax_obs(parallaxes_obs[i],pTrue[i],sigma_parallaxes_obs[i])
		holder+=lnpriorparallax_true(pTrue[i])
		holder+=lnpriorMagnitude_abs_true(MTrue[i],muMTrue,varMTrue**0.5)
	holder+=lnpriormu_M(muMTrue)
	holder+=lnpriorsigma_M(varMTrue**0.5)
	return holder

#lnposterior([10 for i in range(802)],[5 for i in range(802)],

#Setting up the MCMC method, starting position will be an educated guess from the observations we have.
ndim = 2*N+2
nwalkers = ndim * 5

#Some negative parallaxes cause NaN. Replace the values with -15. (-15 was smallest number I found in the list)
initialMTrue = magnitudes_app_obs + 5*np.log10(parallaxes_obs) + 5
for i in range(len(initialMTrue)):
	if np.isnan(initialMTrue[i]):
		initialMTrue[i]=-15

initialpTrue = parallaxes_obs
initialmuMTrue = np.mean(initialMTrue)
initialvarMTrue = np.var(initialMTrue)

#parameter order is: [mu_MTrue,var_MTrue,parallaxes_true[0],...,parallaxes_true[399],Mag_true[0],...,Mag_true[399]]
initialParameters = np.concatenate((np.array([initialmuMTrue, initialvarMTrue]),initialpTrue, initialMTrue))

#Initialize the walkers around this initialguess, but each in a slightly different position. This does deviate from the initialization position that Anthony uses.
initialPos = [initialParameters + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

#args are the observations/Data that is known. The parameters are the positions of the walkers. The order of the inputvariables of the lnposterior needs to be: parameters, observations. (which it is)
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnposterior, args=(magnitudes_app_obs,parallaxes_obs,sigma_magnitudes_app_obs,sigma_parallaxes_obs))

#We make 250k "burn-in" steps, after which we store the position, prob and state (state is the state of the randomnumber generator). Reset and restart from "burn-in" position using the left off state of the random number generator. Store only every 150th sample
pos,prob,state = sampler.run_mcmc(initialPos, 250000)
sampler.reset()
sampler.run_mcmc(pos,10**6,rstate0=state, thin=150)

#Flatchain makes it a 2d array where the rows are the runs and the columns the parameters/dimensions in that run
chain = sampler.flatchain

#Some first analysis using this chain.
meanAbsoluteMagnitudeSamples = chain[:,0].flatten()
varAbsoluteMagnitudeSamples = chain[:,1].flatten()
estimatedAbsMag=meanAbsoluteMagnitudeSamples.mean()
errorEstimatedAbsMag=meanAbsoluteMagnitudeSamples.std()
estimatedVarMag=varAbsoluteMagnitudeSamples.mean()
errorEstimatedVarMag=varAbsoluteMagnitudeSamples.std()
print("emcee estimates")
print("mu_M={:4.2f}".format(estimatedAbsMag)+" +/- {:4.2f}".format(errorEstimatedAbsMag))
print("sigma^2_M={:4.2f}".format(estimatedVarMag)+" +/- {:4.2f}".format(errorEstimatedVarMag))

#Now we will make a nice corner plot
fig = corner.corner([chain[:,0],chain[:,1]], labels=["$\mu_M$", "$\sigma_M^2$"], truths=[mu_Mag_true, var_Mag_true])
fig.savefig("Ulliemam.png")
plt.gcf().clear()
