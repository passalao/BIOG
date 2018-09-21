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
Obs = netCDF4.Dataset('../../SourceData/Vent/AverageWind.nc')
Wind = Obs.variables['Wind']

# Import Accu data
Obs = netCDF4.Dataset('../../SourceData/Accu/Arthern2006_align.nc', nodata="--")
Accu = Obs.variables['Band1']

# Import Emissivity data
Obs = netCDF4.Dataset('../../SourceData/WorkingFiles/Emissivity_FromMatzler_GaussLegendre.nc')#radientDescent_Scipy4QGIS.nc')
E = Obs.variables['Emissivity']

#Resize the Wind and Accu files
Wind=Wind[:,3:-3]
Accu=Accu[:,1:-1]

'''fig, ax = plt.subplots(nrows=1, ncols=1)
#norm = mpl.colors.Normalize(vmin=10, vmax=30)
norm = mpl.colors.Normalize(vmin=0.95, vmax=1)
cmap = mpl.cm.spectral
print(np.shape(E))
myplot = ax.pcolormesh(E[:,112:224], cmap=cmap, norm=norm)
#cbar = fig.colorbar(myplot, ticks=np.arange(0.90, 1.01, 0.02))
#cbar.set_label('Emissivity', rotation=270)
#cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
plt.show()'''


#Selection : East Antarctica, or West Antarctica
Zone="East"#"West"
if Zone=="East":
    E=E[:,112:224]
    Wind = Wind[:,112:224]
    Accu = Accu[:,112:224]

if Zone=="West":
    E=E[:,0:112]
    Wind = Wind[:,0:112]
    Accu = Accu[:, 0:112]

#Selection : Accumulation or Wind
Param="Accu" #"Accu" or "Wind"

Wind1D=np.reshape(Wind, (1, np.size(Wind)))
Accu1D=np.reshape(Accu, (1, np.size(Accu)))
Accu1D[Accu1D>1e+38]=0.0

E1D=np.reshape(E, (1, np.size(E)))
Wind1D=Wind1D[E1D<=1]
Accu1D=Accu1D[E1D<=1]
E1D=E1D[E1D<1]

Wind1D = Wind1D[E1D > 0.95]
Accu1D = Accu1D[E1D > 0.95]
E1D=E1D[E1D>0.95]

if Param=="Wind":
    x=Wind1D
    minx=3
    maxx=12
if Param=="Accu":
    x=Accu1D
    minx=0
    maxx=0.3

y=E1D
y=y[x>0]
x=x[x>0]

'''cmap='jet'
norm = mpl.colors.Normalize(vmin=0.95,vmax=1.)
ax=plt.scatter(Wind1D, Accu1D, s=1.5, c=E1D, cmap=cmap, norm=norm)
plt.colorbar(ax)
plt.xlim(10,28)
plt.xlabel('Wind (m/s)')
plt.ylabel('Accu (m/a)')
plt.ylim(0,0.5)
plt.plot()
plt.show()'''


##############################################################
#Plot kernel density
import pandas as pd
import seaborn as sns
sns.set(style="white")

# Show the joint distribution using kernel density estimation
x1 = pd.Series(x)#, name="Wind speed (m/s)")
x2 = pd.Series(y)#, name="Emissivity")
g = sns.kdeplot(x1, x2, shade=True)# ,cmap=Color)

#Plot the regression line, on a selection of points
if Param=="Wind":
    '''y = y[x < 18]
    x = x[x < 18]
    y = y[x > 10]
    x = x[x > 10]'''
    y = y[x < 10]
    x = x[x < 10]
    y = y[x > 5]
    x = x[x > 5]
if Param=="Accu":
    y=y[x<0.18]
    x=x[x<0.18]
    y=y[x>0.1]
    x=x[x>0.1]

from sklearn import linear_model
regr = linear_model.RANSACRegressor()
regr.fit(x[:,np.newaxis], y)
#print("Coef : ", regr.estimator_.coef_)
x_test = np.linspace(10, 24, 100)
x_test = np.linspace(0, 12, 100)
#x_test = np.linspace(0, 0.5, 100)

plt.plot(x_test, regr.predict(x_test[:,np.newaxis]), color='blue', linewidth=1)

#Design of the plot features
plt.xlim(minx, maxx)
plt.ylim(0.95,1)
plt.xlabel('Wind speed (m/s)', fontsize=17)
plt.xlabel('Accumulation (m/a)', fontsize=17)
plt.ylabel('Emissivity', fontsize=17)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.text(11,.955, Zone, size=17)
plt.show()