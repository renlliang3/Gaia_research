import numpy as np
from scipy import stats
from mayavi import mlab

mu=np.array([1,10,20])
# Let's change this so that the points won't all lie in a plane...
sigma=np.matrix([[20,10,10],
                 [10,25,1],
                 [10,1,50]])

data=np.random.multivariate_normal(mu,sigma,1000)
values = data.T

kde = stats.gaussian_kde(values)

# Create a regular 3D grid with 50 points in each dimension
xmin, ymin, zmin = data.min(axis=0)
xmax, ymax, zmax = data.max(axis=0)
xi, yi, zi = np.mgrid[xmin:xmax:50j, ymin:ymax:50j, zmin:zmax:50j]

# Evaluate the KDE on a regular grid...
coords = np.vstack([item.ravel() for item in [xi, yi, zi]])
density = kde(coords).reshape(xi.shape)

# Visualize the density estimate as isosurfaces
mlab.contour3d(xi, yi, zi, density, opacity=0.5)
mlab.axes()
mlab.show()
