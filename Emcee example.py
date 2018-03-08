import numpy as np
import matplotlib.pyplot as plt
import emcee

def lnprob(x, mu, icov):
    diff = x-mu
    return -np.dot(diff,np.dot(icov,diff))/2.0

ndim = 50

means = np.random.rand(ndim) #Compute 50 random means

cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim)) #Creating random entries between -0.5 and +0.5 in ndim square matrix
cov = np.triu(cov) #leaves upper triangle of the matrix as is and fills the bottom triangle with zeros
cov += cov.T - np.diag(cov.diagonal()) #we add the tranpose of the matrix so that it is symmetric and we replace diag elements with zero
cov = np.dot(cov,cov) #replace matrix by dot product of the matrix with itself for symmetric matrix with positive diag elements

icov = np.linalg.inv(cov) #inverse of the cov matrix for ez calculations later on

nwalkers = 250 #we will use 250 walkers for our task
p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim)) #250x50 matrix with randomly inialitized starting position for each walker

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[means, icov]) #args makes us have to call like this: lnprob(p,means,icov)

#First we do 100 burn-in steps and we save the values for pos, prob and state after these 100 steps
pos, prob, state = sampler.run_mcmc(p0, 100)
sampler.reset() #We reset so that we can start fresh later with the pos we saved

sampler.run_mcmc(pos, 1000) #Production run of 1000 steps

for i in range(ndim):
    plt.figure()
    plt.hist(sampler.flatchain[:,i], 100, color="k", histtype="step") #sampler.flatchain are the 250k positions in 50-dim space
    plt.title("Dimension {0:d}".format(i))

plt.show()

print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
