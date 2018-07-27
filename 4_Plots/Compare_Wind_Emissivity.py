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
Obs = netCDF4.Dataset('../../SourceData/Vent/MaxWind.nc')
Wind = Obs.variables['Wind']

# Import Emissivity data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/WorkingFiles/Emissivity_FromMatzler.nc')#radientDescent_Scipy4QGIS.nc')
E = Obs.variables['Emissivity']

#Resize the Wind file
Wind=Wind[:,3:-3]
print(np.shape(Wind))
Wind1D=np.reshape(Wind, (1, np.size(Wind)))
E1D=np.reshape(E, (1, np.size(E)))

Wind1D2=Wind1D[E1D<1]
E1D2=E1D[E1D<1]

x=Wind1D2[E1D2>0.9]
y=E1D2[E1D2>0.9]

import pandas as pd
import seaborn as sns
sns.set(style="white")

# Generate a random correlated bivariate dataset
rs = np.random.RandomState(5)
mean = [0, 0]
cov = [(1, .5), (.5, 1)]
#x1, x2 = rs.multivariate_normal(mean, cov, 500).T

#Linear regression
#y=E1D3[Wind1D3<16]
#x=Wind1D3[Wind1D3<16]
#y=y[x>10]
#x=x[x>10]

#Normalize the data
Normx=(x-min(x))/(max(x)-min(x))
Normy=(y-min(y))/(max(y)-min(y))


from sklearn import linear_model
regr = linear_model.LinearRegression()
regr.fit(Normx[:,np.newaxis], Normy)
print(regr.coef_)

x_test = np.linspace(np.min(Normx), np.max(Normx), 100)
plt.plot(x_test, regr.predict(x_test[:,np.newaxis]), color='blue', linewidth=3)
#plt.show()
plt.scatter(Normx,Normy, s=1e-1)
plt.show()

'''x1 = pd.Series(x, name="Wind speed (m/s)")
x2 = pd.Series(y, name="Emissivity")

# Show the joint distribution using kernel density estimation
g = sns.kdeplot(x1, x2, shade=True)
plt.xlim(10,24)
#plt.xlim(3,10)
plt.ylim(0.95,1)
plt.savefig("../../Article/img/WindScatter.eps")
plt.show()'''

'''plt.scatter(x1,x2,s=1)
plt.xlim(10,24)
#plt.xlim(3,10)
plt.ylim(0.95,1)
plt.show()'''