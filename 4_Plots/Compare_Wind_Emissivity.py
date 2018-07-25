from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import NC_Resources as ncr
import scipy.optimize as opt
#np.set_printoptions(threshold=np.nan)


# Import Wind data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/Vent/AveragedWind.nc')
Wind = Obs.variables['Wind']

# Import Emissivity data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/WorkingFiles/Emissivity_FromGradientDescent_Scipy4QGIS.nc')
E = Obs.variables['Emissivity']

#Resize the Wind file
Wind=Wind[:,3:-3]
print(np.shape(Wind))
Wind1D=np.reshape(Wind, (1, np.size(Wind)))
E1D=np.reshape(E, (1, np.size(E)))

Wind1D2=Wind1D[E1D<1]
E1D2=E1D[E1D<1]

Wind1D3=Wind1D2[E1D2>0.9]
E1D3=E1D2[E1D2>0.9]

import pandas as pd
import seaborn as sns
sns.set(style="white")

# Generate a random correlated bivariate dataset
rs = np.random.RandomState(5)
mean = [0, 0]
cov = [(1, .5), (.5, 1)]
x1, x2 = rs.multivariate_normal(mean, cov, 500).T
x1 = pd.Series(Wind1D3, name="Wind speed (m/s)")
x2 = pd.Series(E1D3, name="Emissivity")

# Show the joint distribution using kernel density estimation
g = sns.kdeplot(x1, x2, shade=True)
plt.xlim(3,12)
plt.ylim(0.95,1)
plt.savefig("../../Article/img/WindScatter.eps")
plt.show()