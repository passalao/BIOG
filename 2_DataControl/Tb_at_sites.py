#!/usr/bin/python
# -*- coding: cp1252 -*-

###############################################################
#Plot the temperature profiles at drill sites                 #
###############################################################
from netCDF4 import Dataset
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import pyproj
import NC_Resources as ncr
from numpy import loadtxt

##############################################################
#Import data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_57.5deg_xy.nc')
Tb_Obs = Obs.variables['BT_V']
X=Obs.variables['x_ease2']
Y=Obs.variables['y_ease2']
Lx=25e3
Ly=25e3

#Import GRISLI data
#GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4Ts.nc')
GRISLI = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/Corrected_Tz_MappedonSMOS.nc')
H = GRISLI.variables['H']
Zeta = GRISLI.variables['Zeta']
'''Zeta=np.zeros((201,225,21))
for k in np.arange(0,21,1):
    Zeta[:,:,k]=1-0.05*k'''

#Zeta=np.linspace(1,0,0.05)
#Zeta=Zeta[:,100,100]
#print(Zeta)
TzGR = GRISLI.variables['T']
nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(GRISLI)

#Drill sites
Sites=['DomeC', 'Vostok', 'DomeFuji', 'EDML','Byrd','LawDome']#, 'Berkner']#, 'SipleDome',]
Lon=[123.3952,106.7114,39.7222,0.05161,-119.31,112.8067]#, -45.6783]#, -149.2426,]
Lat=[-75.1017,-78.4719,-77.3088,-75.00,-80.01,-66.7391]#, -79.5483]#, -81.6620,]
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj(init="EPSG:6932")
Xs, Ys = pyproj.transform(wgs84, StereoPol, Lon, Lat)
Xpix=(Xs-X[0][0])//Lx+1
Ypix=(Ys-Y[-1][-1])//Ly+1

##############################################################
#Plot
Ts=np.zeros(np.shape(Sites))
Tb=np.zeros(np.shape(Sites))

k=0
TsAtSites=[]
colors=['r','b','g','k','orange','purple','pink',"cyan","gray"]
for site,c in zip(Sites, colors):
    j=int(Ypix[k])
    i=int(Xpix[k])
    Data = loadtxt("../../SourceData/Temperatures/"+str(site)+".csv", comments="#", delimiter=",",unpack=False)
    Tz_Obs = Data[:, 1]
    depth = Data[:, 0]
    Ts[k] = Tz_Obs[0]
    Tb[k]=Tb_Obs[0, j, i]

    plt.plot(Tz_Obs, depth, label=site, color=c, linewidth=1)
    [plt.scatter(Tb_Obs[0,j,i]-273.15, 0,color=c) for i,j,c in zip(Xpix, Ypix, colors)]
    TsAtSites.append(Tz_Obs[0])

    # Compute the error of GRISLI
    print(j,i)
    depthini=(1 - Zeta[j,i]) * H[j,i]
    TzGR_interp=np.interp(depth, depthini, TzGR[j,i])
    print("Mean error of GRISLI at", site, ':',(np.dot(TzGR_interp-Tz_Obs[:],TzGR_interp-Tz_Obs[:])/np.size(TzGR_interp))**0.5, "K")

    k=k+1

#Avec correction pour Ts:
#[plt.plot(TzGR[j,i]-TzGR[j,i][0]+TsAtSites[k], (1 - Zeta[j,i]) * H[j,i], '--',color=c, linewidth=1) for i, j, k, site, c in zip(Xpix,Ypix, np.arange(0,np.size(Xpix),1),Sites, colors)]
#Sans correction pour Ts :
[plt.plot(TzGR[j,i], (1 - Zeta[int(j),int(i)]) * H[int(j),int(i)], '--',color=c, linewidth=1) for i, j, k, site, c in zip(Xpix,Ypix, np.arange(0,np.size(Xpix),1),Sites, colors)]

plt.legend()
plt.xlabel("Ice temperature (°C)", fontsize=13)
plt.ylabel("Depth (m)",  fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.gca().invert_yaxis()
plt.xlim(-60,0)
plt.legend()
plt.grid()
plt.show()