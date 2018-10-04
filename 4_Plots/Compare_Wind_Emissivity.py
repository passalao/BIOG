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
#np.set_printoptions(threshold=np.nan)

##############################################################
#Data import

# Import Wind data
Obs = netCDF4.Dataset('../../SourceData/Vent/AverageWind.nc')
Wind = Obs.variables['Wind']

# Import Accu data
Obs = netCDF4.Dataset('../../SourceData/Accu/Arthern2006_align.nc', nodata="--")
Accu = Obs.variables['Band1']

# Import Emissivity data
Obs = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/Emissivity_FromMatzler_GaussLegendre.nc')#radientDescent_Scipy4QGIS.nc')
E = Obs.variables['Emissivity']
TETs = Obs.variables['TE-Ts']
TE = Obs.variables['Teff']

#Import TbV data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg.nc')
TbV = Obs.variables['BT_V']

#Import TbH data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbH_52.5deg.nc')
TbH = Obs.variables['BT_H']

#Resize the Wind and Accu files
Wind=Wind[:,3:-3]
Accu=Accu[:,1:-1]
Accu=Accu[::-1,:]

'''fig, ax = plt.subplots(nrows=1, ncols=1)
#norm = mpl.colors.Normalize(vmin=3, vmax=10)
#norm = mpl.colors.Normalize(vmin=0.95, vmax=1)
#norm = mpl.colors.Normalize(vmin=0, vmax=0.05)

cmap = mpl.cm.spectral
print(np.shape(E))
myplot = ax.pcolormesh(E, cmap=cmap, norm=norm)
#cbar = fig.colorbar(myplot, ticks=np.arange(0.90, 1.01, 0.02))
#cbar.set_label('Emissivity', rotation=270)
#cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
plt.show()'''


#Selection : East Antarctica, or West Antarctica
Zone="East"#"West"
'''if Zone=="East":
    E=E[:,112:224]
    TETs=TETs[:,112:224]
    Wind = Wind[:,112:224]
    Accu = Accu[:,112:224]

if Zone=="West":
    E=E[:,0:112]
    TETs=TETs[:,0:112]
    Wind = Wind[:,0:112]
    Accu = Accu[:, 0:112]'''

#Selection : Accumulation or Wind
Param="Wind" #"Accu" or "Wind"

Wind1D=np.reshape(Wind, (1, np.size(Wind)))
Accu1D=np.reshape(Accu, (1, np.size(Accu)))
TETs1D=np.reshape(TETs, (1, np.size(Accu)))
TE1D = np.reshape(TE, (1, np.size(Accu)))
TbH1D = np.reshape(TbH, (1, np.size(Accu)))
TbV1D = np.reshape(TbV, (1, np.size(Accu)))
Accu1D[Accu1D>1e+38]=0.0


E1D=np.reshape(E, (1, np.size(E)))
E1D[E1D>1]=1

'''Wind1D=Wind1D[E1D<=1]
Accu1D=Accu1D[E1D<=1]
TETs1D=TETs1D[E1D<=1]
E1D=E1D[E1D<1]'''

'''Wind1D = Wind1D[E1D > 0.95]
Accu1D = Accu1D[E1D > 0.95]
TETs1D=TETs1D[E1D>0.95]
E1D=E1D[E1D>0.95]'''

'''Wind1D = Wind1D[TETs1D > 0]
Accu1D = Accu1D[TETs1D > 0]
E1D=E1D[TETs1D > 0]
TETs1D=TETs1D[TETs1D > 0]'''

if Param=="Wind":
    minx=3
    maxx=12
if Param=="Accu":
    minx=0
    maxx=0.3

#y=y[x>0]
#x=x[x>0]

##############################################################
#Plot kernel density
import pandas as pd
import seaborn as sns
sns.set(style="white")

# Show the joint distribution using kernel density estimation
#x1 = pd.Series(x)#, name="Wind speed (m/s)")
#x2 = pd.Series(y)#, name="Emissivity")
#g = sns.kdeplot(x1, x2, shade=True, alpha=0.5)# ,cmap=Color)

#Cleaning of dummy values
#Wind1D=Wind1D[Accu1D>0]
#E1D=E1D[Accu1D>0]
#Accu1D=Accu1D[Accu1D>0]

'''Wind1D=Wind1D[E1D>0.9]
Accu1D=Accu1D[E1D>0.9]
TE1D=TE1D[E1D>0.9]
E1D=E1D[E1D>0.9]'''

'''y=E1D
x = Wind1D

#Scatterplot
cmap='coolwarm'
norm = mpl.colors.Normalize(vmin=0, vmax=0.3)
ax=plt.scatter(x, y, c=Accu1D, s=5, linewidth=0, alpha=0.7, norm=norm, cmap=cmap)
cbar=plt.colorbar(ax)
cbar.set_label('Accumulation rate (m/a)', rotation=90,fontsize=14, labelpad=10)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation='vertical')
'''
'''Accu1D=Accu1D[x<1]
y=y[x<1]
x=x[x<1]'''

'''#Scatterplot
cmap='rainbow'
norm = mpl.colors.Normalize(vmin=0.03, vmax=0.3)
ax=plt.scatter(E1D, TbV1D-TbH1D, c=Accu1D, s=5, linewidth=0, alpha=0.7, norm=norm, cmap=cmap)
cbar=plt.colorbar(ax)
cbar.set_label('Accumulation rate (m/a)', rotation=90,fontsize=14, labelpad=10)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation='vertical')
x=x[TETs1D > 0]
y=y[TETs1D > 0]
plt.xlim(0.93,1.00)
plt.ylim(20,60)
plt.show()'''

'''#Plot the regression line, on a selection of points
if Param=="Wind":
    y = y[x < 9]
    x = x[x < 9]
    y = y[x > 4]
    x = x[x > 4]
if Param=="Accu":
    y=y[x<0.10]
    x=x[x<0.10]
    #y=y[x>0.12]
    #x=x[x>0.12]'''

'''from sklearn import linear_model
y=E1D[Accu1D>0.2]
x=Wind1D[Accu1D>0.2]
regr = linear_model.RANSACRegressor()
regr.fit(x[:,np.newaxis], y)
print("Coef : ", regr.estimator_.coef_)
x_test = np.linspace(0, 12, 100)
#x_test = np.linspace(0, 0.5, 100)
plt.plot(x_test, regr.predict(x_test[:,np.newaxis]), color='red', linewidth=1)

y=E1D[Accu1D<0.10]
x=Wind1D[Accu1D<0.10]
regr = linear_model.RANSACRegressor()
regr.fit(x[:,np.newaxis], y)
print("Coef : ", regr.estimator_.coef_)
#x_test = np.linspace(10, 24, 100)
x_test = np.linspace(0, 12, 100)
#x_test = np.linspace(0, 0.5, 100)
plt.plot(x_test, regr.predict(x_test[:,np.newaxis]), color='blue', linewidth=1)


y=E1D[Accu1D<0.2]
x=Wind1D[Accu1D<0.2]
Accu1D=Accu1D[Accu1D<0.2]
y=y[Accu1D>0.1]
x=x[Accu1D>0.1]
regr = linear_model.RANSACRegressor()
regr.fit(x[:,np.newaxis], y)
print("Coef : ", regr.estimator_.coef_)
x_test = np.linspace(0, 12, 100)
#x_test = np.linspace(0, 0.5, 100)
plt.plot(x_test, regr.predict(x_test[:,np.newaxis]), color='gray', linewidth=1)'''


'''#Design of the plot features
plt.xlim(minx, maxx)
plt.ylim(0.92,1.00)
#plt.xlabel(Param+" "+unit, fontsize=14)
plt.xlabel("Wind speed (m/s)", fontsize=14)
plt.ylabel('Emissivity', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.text(4,.955, Zone, size=17)
#plt.show()'''

#Here compute the correlation between emissivity and TE
T=np.linspace(-60,-5,12)

for t in T:
    print(t)
    min=t
    max=t+5
    x=E1D[TE1D<273.15+max]
    y=TE1D[TE1D<273.15+max]

    x=x[y>273.15+min]
    y=y[y>273.15+min]

    #x=E1D#[0,:]
    #y=TE1D#[0,:]
    '''y=y[x>0]#[0,:]
    x=x[x>0]#[0,:]
    y=y[x<1]#[0,:]
    x=x[x<1]'''

    '''cmap='rainbow'
    norm = mpl.colors.Normalize(vmin=0.03, vmax=0.3)
    ax=plt.scatter(x,y, s=1, alpha=0.5)#, norm=norm, cmap=cmap)
    #cbar=plt.colorbar(ax)
    plt.show()'''

    m = np.stack((x,y), axis=0)
    sigmax = np.std(x)
    sigmay = np.std(y)
    #print(min(x), max(x), min(y), max(y))
    #print(sigmax, sigmay)

    #print((sum((x-np.mean(x))*(y-np.mean(y)))/sigmax/sigmay/np.size(x))**0.5)
    print((abs(np.cov(m)[0, 1]) / sigmax / sigmay)**0.5)
    #print((np.cov(m)[1, 1] / sigmay / sigmay)**0.5)