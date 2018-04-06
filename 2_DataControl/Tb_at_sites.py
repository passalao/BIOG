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

Sites=['DomeC', 'Vostok', 'DomeFuji', 'EDML','SipleDome','Byrd','LawDome']#, 'T-UC', 'T-AIS']
Lon=[123.3952,106.7114,39.7222,0.05161,-149.2426,-119.31,112.8067]#,-138.372,-138.946]
Lat=[-75.1017,-78.4719,-77.3088,-75.00,-81.6620,-80.01,-66.7391]#,-83.679,-83.4619]
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
Tb=np.zeros(np.shape(Sites))
DeltaTTiuri=np.zeros(np.shape(Sites))
DeltaTMatzler=np.zeros(np.shape(Sites))
DeltaTFirn=np.zeros(np.shape(Sites))

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
    #Tb_modobs=[t*0.98 for t in Tb_modobs]
    Error[k]=Tb_modobs-Tb_Obs[0,j,i]
    Height[k]=depth[-1]
    Ts[k]=Tz_Obs[0]
    Tb[k]=Tb_Obs[0, j, i]

    # Compute T at d=466 and d=1750
    temp=np.interp([0, 100, 200, 300, 400, 466, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1750,
         1800], depth, Tz_Obs)
    DeltaTTiuri[k]=temp[5]-Ts[k]
    DeltaTMatzler[k] = temp[19] - Ts[k]
    DeltaTFirn[k] = temp[1] - Ts[k]
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

'''plt.scatter(Error[:,0],Error[:,1], c=Height, cmap=cmap)
[plt.text(et-1,em-1,site) for et,em,site in zip(Error[:,0], Error[:,1],Sites)]
cbar=plt.colorbar()
cbar.set_label('Ice thickness (m)', rotation=270, labelpad=15)
plt.xlabel("Error with Tiuri")
plt.ylabel("Error with Mätzler")
plt.xlim(-2,18)
plt.ylim(-2,18)'''
'''#plt.scatter(Error[:,0],DeltaTFirn, c=Height, cmap=cmap)
#[plt.text(et-1,dt-0.1,site) for et,dt,site in zip(Error[:,0], DeltaTFirn,Sites)]
#plt.scatter(Error[:,1],DeltaTMatzler, c=Height, cmap=cmap)
#[plt.text(et-1,dt-0.5,site) for et,dt,site in zip(Error[:,1], DeltaTMatzler,Sites)]
#plt.scatter(Error[:,0],DeltaTTiuri, c=Height, cmap=cmap)
#[plt.text(et-1,dt-0.5,site) for et,dt,site in zip(Error[:,0], DeltaTTiuri,Sites)]'''
plt.scatter(Ts+273.15,Tb, c=Height, cmap=cmap)
[plt.text(ts-1+273.15,tb-3,site) for ts,tb,site in zip(Ts,Tb,Sites)]
cbar=plt.colorbar()
cbar.set_label('Ice thickness (m)', rotation=270, labelpad=15)
#plt.xlabel("Error with Tiuri")
#plt.ylabel("T(466 m)-Ts")
#plt.xlabel("Error with Mätzler")
#plt.ylabel("T(1750 m)-Ts")
#plt.xlabel("Error with Tiuri")
#plt.ylabel("T(100 m)-Ts")
plt.xlabel("Ts (K)")
plt.ylabel("Tb SMOS (K)")

#plt.xlim(-2,15)
#plt.ylim(-2,15)
plt.grid(which='both')
plt.axis("equal")
#plt.autoscale(False)
#plt.gca().autoscale_view()
plt.xlim(210,260)
plt.ylim(210,260)
plt.plot([200,270],[200,270],'-', c='r', lw=0.5)
#plt.plot([-5,18],[-5,18],'-', c='r', lw=0.5)
'''from sklearn import linear_model
regr = linear_model.LinearRegression()
regr.fit(Ts[:,np.newaxis], Tb)
Ts_test = np.linspace(np.min(Ts), np.max(Ts), 100)
plt.plot(Ts_test+273.15, regr.predict(Ts_test[:,np.newaxis]), color='blue', linewidth=1)'''
plt.show()

'''#Normalization
Error[:,1]=(Error[:,1]-min(Error[:,1]))/(max(Error[:,1])-min(Error[:,1]))
Error[:,0]=(Error[:,0]-min(Error[:,0]))/(max(Error[:,0])-min(Error[:,0]))'''

cmap = plt.get_cmap('inferno')
colors = [cmap(i) for i in Error[:,1]]

#Plot
i=0
for site in Sites:
    #color=colors[i]
    Data = loadtxt("../../SourceData/Temperatures/"+str(site)+".csv", comments="#", delimiter=",",unpack=False)
    Tz_Obs = Data[:, 1]
    depth = Data[:, 0]
    plt.plot(Tz_Obs, depth, label=site)# color=color)
    plt.legend()
    i=i+1

#cbar = plt.colorbar()
#cbar.set_label('Error (K)', rotation=270, labelpad=15)
plt.xlabel("Ice temperature (K)")
plt.ylabel("Depth (m)")
plt.gca().invert_yaxis()
plt.grid()
plt.show()

