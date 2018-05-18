from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

# Import temperature data
Emiss = netCDF4.Dataset('../../SourceData/WorkingFiles/Emissivity.nc')
E = np.array(Emiss.variables['Emissivity'])

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']
Tb=Tb[0]
Lx=25e3
Ly=25e3

#Sites
SitesS=['DomeC', 'Vostok', 'DomeFuji', 'EDML','SipleDome','Byrd','LawDome']#, 'T-UC', 'T-AIS']
Lon=[123.3952,106.7114,39.7222,0.05161,-149.2426,-119.31,112.8067]#,-138.372,-138.946]
Lat=[-75.1017,-78.4719,-77.3088,-75.00,-81.6620,-80.01,-66.7391]#,-83.679,-83.4619]
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPolS = pyproj.Proj(init="EPSG:6932")
XsS, YsS = pyproj.transform(wgs84, StereoPolS, Lon, Lat)
XpixS=(XsS-X[0][0])//Lx+1
YpixS=(YsS-Y[-1][-1])//Ly+1

#Emissivity computed with Tiuri
ETiuri=[0.985,0.995,0.98,0.992,0.964,0.979,0.987, 0.988]
EGRISLI=[E[int(j),int(i)] for i,j in zip(XpixS, 201-YpixS)]

[plt.scatter(et, eg, label=site) for site, et, eg in zip(SitesS, ETiuri, EGRISLI)]
plt.plot([0,1],[0,1], c='r')
plt.axis("equal")
plt.xlim(0.96,1)
plt.ylim(0.96,1)
plt.xlabel("Emissivity from boreholes and Tiuri")
plt.ylabel("Emissivity from GRISLI and regression model")

plt.legend()
plt.grid()
plt.show()
