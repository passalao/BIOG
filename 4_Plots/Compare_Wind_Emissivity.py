###############################################################
#Plot the joint distribution of emissivity and wind speed     #
###############################################################
from netCDF4 import Dataset
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import NC_Resources as ncr
import scipy.optimize as opt

##############################################################
#Data import

# Import Wind data
Obs = netCDF4.Dataset('../../SourceData/Vent/MaxWind.nc')
Wind = Obs.variables['Wind']

# Import Emissivity data
Obs = netCDF4.Dataset('../../SourceData/WorkingFiles/Emissivity_FromMatzler.nc')#radientDescent_Scipy4QGIS.nc')
E = Obs.variables['Emissivity']

#Resize the Wind file
Wind=Wind[:,3:-3]
Wind1D=np.reshape(Wind, (1, np.size(Wind)))
E1D=np.reshape(E, (1, np.size(E)))
Wind1D=Wind1D[E1D<1]
E1D=E1D[E1D<1]
x=Wind1D[E1D>0.95]
y=E1D[E1D>0.95]

##############################################################
#Plot kernel density
import pandas as pd
import seaborn as sns
sns.set(style="white")

x1 = pd.Series(x, name="Wind speed (m/s)")
x2 = pd.Series(y, name="Emissivity")

# Show the joint distribution using kernel density estimation
g = sns.kdeplot(x1, x2, shade=True)
plt.xlim(10,24)
plt.ylim(0.95,1)

#plot the regression line, on a selection of points
y=y[x<16]
x=x[x<16]
from sklearn import linear_model
regr = linear_model.RANSACRegressor()
regr.fit(x[:,np.newaxis], y)
x_test = np.linspace(np.min(x), np.max(x), 100)
plt.plot(x_test, regr.predict(x_test[:,np.newaxis]), color='blue', linewidth=1)
#inliers and outliers of the RANSAC process
#inlier_mask=regr.inlier_mask_
#outlier_mask = np.logical_not(inlier_mask)
#plt.scatter(x[outlier_mask], y[outlier_mask], color='gold', marker='.',label='Outliers')

plt.show()