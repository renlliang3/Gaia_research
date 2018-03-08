import matplotlib.pyplot as plt
import numpy as np

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


#parallax=gaussian(np.linspace(0,10,1000),5,2)

def plotdistribution(mean,sd):
	parallax=np.random.normal(mean,sd,1000)
	distance=1/parallax
	plt.subplot(1, 2, 1)
	plt.hist(parallax)
	plt.subplot(1, 2, 2)
	plt.hist(distance)
	plt.savefig('mean='+str(mean)+' sd='+str(sd))

def distribution(mean,sd):
	parallax=np.random.normal(mean,sd,1000)
	distance=1/parallax
	difference=abs(np.mean(distance)-1/np.mean(parallax))
	return difference

def drawerror(k):
	x=10*np.array(range(k))+10
	y=[distribution(200,i) for i in x]
	plt.plot(x,y)
	plt.xlabel('Standard deviation')
	plt.ylabel(r'$|E[1/\omega_0]-1/\omega|$')
	plt.savefig('Error growing')

plotdistribution(100,90)
