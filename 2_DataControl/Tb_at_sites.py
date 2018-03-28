#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
#import mpl_toolkits.basemap.pyproj as pyproj
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
from numpy import loadtxt

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
Lon = Obs.variables['lon']
Lat = Obs.variables['lat']
X=Obs.variables['x_ease2']
Y=Obs.variables['y_ease2']
Mask = Obs.variables['mask']
#Pixel size
Lx=25e3
Ly=25e3

#Import GRISLI data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
H = GRISLI.variables['H']
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']

Sites=['DomeC', 'Vostok', 'DomeFuji', 'SipleDome','Byrd','LawDome', 'T-UC', 'T-AIS']
Lon=[123.3952,106.7114,39.7222,-149.2426,-119.31,112.8067,-138.372,-138.946]
Lat=[-75.1017,-78.4719,-77.3088,-81.6620,-80.01,-66.7391,-83.679,-83.4619]
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj(init="EPSG:6932")  #
Xs, Ys = pyproj.transform(wgs84, StereoPol, Lon, Lat)
Xpix=(Xs-X[0][0])//Lx+1
Ypix=(Ys-Y[-1][-1])//Ly+1

nbfields=7
OutData=np.zeros((np.shape(Sites)[0],nbfields))

#f=open("../../OutputData/RT_Sites_"+str(BIOG.var.Perm)+'_'+str(BIOG.var.NbLayers)+"l.csv",'a')
f=open("../../OutputData/RT_Sites.csv",'a')
f.write('\n')
f.write('############################# \n')
f.write('Number_of_vertical_layers: '+str(BIOG.var.NbLayers)+"\n")
f.write('Permittivity_model: '+str(BIOG.var.Perm)+"\n")
f.write('Number_of_streams: '+str(BIOG.var.NbStreams)+"\n")
f.write('\n')
f.write('Site Lon Lat Thickness TbSMOS Ts_LocalObs Tb_DMRT-ML Error_Obs\n')

Error=np.zeros((np.shape(Sites)[0],2))
Height=np.zeros(np.shape(Sites))
Ts=np.zeros(np.shape(Sites))

Perm=['Tiuri', 'Matzler']
k=0
for site in Sites:
    print(site)
    j=Ypix[k]
    i=Xpix[k]
    Tz_gr_at_Point = Tz_gr[j,i, :]
    Data = loadtxt("../../SourceData/Temperatures/"+str(site)+".csv", comments="#", delimiter=",",unpack=False)
    Tz_Obs = Data[:, 1]
    depth = Data[:, 0]

    # Compute the density for the whole profile
    Tb_modobs = [BIOG.fun.GetTb(Tz_Obs[:], depth[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams,
                            p, BIOG.var.RTModel,0) for p in Perm]

    #Tb_modgr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams,
    #                       BIOG.var.Perm, BIOG.var.RTModel,0)
    Error[k]=Tb_modobs-Tb_Obs[0,j,i]
    Height[k]=depth[-1]
    Ts[k]=Tz_Obs[0]

    '''OutData[k,0]=Lon[k]
    OutData[k,1] = Lat[k]
    OutData[k,2] = depth[-1]
    OutData[k,3]= Tb_Obs[0,j,i]
    OutData[k,4]=Tz_Obs[0]+273.15
    OutData[k,5] =Tb_modobs
    OutData[k,6] =Error[k]#Tb_modobs-Tb_Obs[0,j,i]

    f.write(str(site)+' ')
    [f.write(str(OutData[k,l])+' ') for l in np.arange(0,nbfields)]
    f.write("\n")'''
    k=k+1


#plt.scatter(Error,Height)

cmap = plt.get_cmap('viridis')
plt.scatter(Error[:,0],Error[:,1], c=Height, cmap=cmap)
[plt.text(et-1,em-1,site) for et,em,site in zip(Error[:,0], Error[:,1],Sites)]
plt.grid()
cbar=plt.colorbar()
cbar.set_label('Ice thickness (m)', rotation=270, labelpad=15)
plt.xlabel("Error with Tiuri")
plt.ylabel("Error with Mätzler")
plt.xlim(-2,18)
plt.ylim(-2,18)
plt.show()

#Normalization
Error=(Error-min(Error))/(max(Error)-min(Error))

cmap = plt.get_cmap('inferno')
colors = [cmap(i) for i in Error]

#Plot
i=0
for site in Sites:
    color=colors[i]
    Data = loadtxt("../../SourceData/Temperatures/"+str(site)+".csv", comments="#", delimiter=",",unpack=False)
    Tz_Obs = Data[:, 1]
    depth = Data[:, 0]
    plt.plot(Tz_Obs, depth, color=color)
    i=i+1

plt.gca().invert_yaxis()
plt.grid()
plt.show()