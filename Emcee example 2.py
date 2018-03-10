import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op
import corner

# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

x_axis=np.linspace(0,10,50)
y_true=x_axis*m_true+b_true

# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
y = m_true*x+b_true
y += np.abs(f_true*y) * np.random.randn(N)
y += yerr * np.random.randn(N)

A = np.vstack((np.ones_like(x), x)).T
C = np.diag(yerr * yerr)
cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))

y_ls=x_axis*m_ls+b_ls

#The likelihood of our Bayes setup model. Which is a Gaussian with something added 
def lnlike(theta, x, y, yerr):
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
m_ml, b_ml, lnf_ml = result["x"]

y_ml=x_axis*m_ml+b_ml

#Here we will start using Bayes theorem to find the probability distribution over all possible parameters

#We will use a uniform prior with a certain parameter range, and zero otherse. This is also the log probability
def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf

#The logarithm of the posterior: add the log of the likelihood the log of the prior
#Remember since we are working with log probability, we have to add probabilities instead of multiplying
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

#Define the starting setup for the MCMC method. Three parameters so ndim=3 and we will use a starting position at the maximum likelihood in parameter space.
ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

#Iniatilize the MCMC method with the starting setup and run it for 500 steps.
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 500)

#It seems that after 50 steps the method has 'burned in' so we will take only the samples after the first 50. Remember the shape of this Chain code:
#which is (walkers, runs of walker, number of dimensions of parameterspace)
#We will reshape this matrix by flattening the first two array dimensions so that we get a 2dim matrix with all the positions in parameterspace
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

#Make a cornerfigure
fig = corner.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],
                      truths=[m_true, b_true, np.log(f_true)])
fig.savefig("triangle.png")
plt.gcf().clear()

xl = np.array([0, 10])
#Plot 100 lines with parameters taken randomly out the samples array
#Plot line for true parameter in red. lw is linewidth. alpha is transparancy of the plotted line
for m, b, lnf in samples[np.random.randint(len(samples), size=100)]:
    plt.plot(xl, m*xl+b, color="k", alpha=0.1)
plt.plot(xl, m_true*xl+b_true, color="r", lw=2, alpha=0.8)
plt.errorbar(x, y, yerr=yerr, fmt=".k")
plt.show()

#because we have the log of we first take the e power for this row
samples[:, 2] = np.exp(samples[:, 2])

#map does something to each value of a list, in the map we define a local function with lambda that will only be used locally in the brackets. v is the input
#parameter of lambda function and after the : is the output value, and lastly you insert what function you want to apply this mapping onto.
#np.percentile computes the median and the bottom and upper part of the variance with 84-16=68. zip rearranges the shape of the array so that the rows
#and columns are switched and like [(param1(50),param1(84),param1(16)),(2(50),2(84),2(16)),(3(50),3(84),3(16))] instead of
#[[1(16),2(16),3(16)],[1(50),2(50),3(50)],[1(84),2(84),3(84)]].
#v[1] will give median, v[2]-v[1] will give the upper difference for the variance, v[1]-v[0] will give the lower difference for the variance.
m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))

print(m_mcmc,b_mcmc,f_mcmc)
