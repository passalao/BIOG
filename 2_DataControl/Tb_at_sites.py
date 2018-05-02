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
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_57.5deg_xy.nc')
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

# Import SMOS data on Greenland
ObsN = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansJJA_TbV_57.5deg_Greenland.nc')
ny_ObsN = ObsN.dimensions['cols'].size
nx_ObsN = ObsN.dimensions['rows'].size
Tb_ObsN = ObsN.variables['BT_V']
XN=ObsN.variables['x_ease2']
YN=ObsN.variables['y_ease2']
#XN=ObsN.variables['x_npstere']
#YN=ObsN.variables['y_npstere']
#Pixel size
LxN=25e3
LyN=25e3

#Import GRISLI data
#GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles_GRISLIini/GRISLIMappedonSMOS.nc')
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_MappedonSMOS.nc')#
H = GRISLI.variables['H']
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Uxbar=GRISLI.variables['UXBAR']
Uybar=GRISLI.variables['UYBAR']
Us_gr=(np.array(Uxbar)**2+np.array(Uybar)**2)**0.5*4/3

SitesS=['DomeC', 'Vostok', 'DomeFuji', 'EDML','SipleDome','Byrd','LawDome']#, 'T-UC', 'T-AIS']
Lon=[123.3952,106.7114,39.7222,0.05161,-149.2426,-119.31,112.8067]#,-138.372,-138.946]
Lat=[-75.1017,-78.4719,-77.3088,-75.00,-81.6620,-80.01,-66.7391]#,-83.679,-83.4619]
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPolS = pyproj.Proj(init="EPSG:6932")
XsS, YsS = pyproj.transform(wgs84, StereoPolS, Lon, Lat)
XpixS=(XsS-X[0][0])//Lx+1
YpixS=(YsS-Y[-1][-1])//Ly+1
print(XpixS, YpixS)

SitesN=['GRIP']
Lon=[-37.62]
Lat=[72.57]
StereoPolN = pyproj.Proj(init="EPSG:6933")
XsN, YsN = pyproj.transform(wgs84, StereoPolN, Lon, Lat)
XsN=[-3629813]
YsN=[7002276]
XpixN=(XsN-XN[0][0])//LxN+1
YpixN=(YsN-YN[0][0])//LyN+1
Sites=SitesS+SitesN
Xpix=np.concatenate([XpixS,XpixN])
Ypix=np.concatenate([YpixS,YpixN])
Emiss=[0.985,1,0.98,0.992,0.964,0.979,0.987, 0.988]

nbfields=7
OutData=np.zeros((np.shape(Sites)[0],nbfields))

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
T500_Obs=np.zeros(np.shape(Sites))
dTdzbar=np.zeros(np.shape(Sites))
DeltaTTiuri=np.zeros(np.shape(Sites))
DeltaTMatzler=np.zeros(np.shape(Sites))
DeltaTFirn=np.zeros(np.shape(Sites))
Perm=['Tiuri', 'Matzler']
k=0
for site in Sites:
    print(site)
    j=Ypix[k]
    i=Xpix[k]
    #Tz_gr_at_Point = Tz_gr[j,i, :]
    Data = loadtxt("../../SourceData/Temperatures/"+str(site)+".csv", comments="#", delimiter=",",unpack=False)
    Tz_Obs = Data[:, 1]
    depth = Data[:, 0]
    #Compute the T gradient on the upper 1000 m
    l=0
    for t,d in zip(Tz_Obs[1:],depth[1:]):
        if d<=500:
            T500_Obs[k]=Tz_Obs[l+1]
        if d <= 1000:
            dTdzbar[k]=dTdzbar[k]+(t-Tz_Obs[l])
        l=l+1

    dTdzbar[k]=dTdzbar[k]/1000

    # Compute the density for the whole profile
    Tb_modobs = [BIOG.fun.GetTb(Tz_Obs[:], depth[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams,
                            p, BIOG.var.RTModel,0) for p in Perm]
    Height[k] = depth[-1]
    Ts[k] = Tz_Obs[0]
    if site!="GRIP":
        Error[k]=np.array(Tb_modobs)*Emiss[k]-Tb_Obs[0,j,i]
        Tb[k]=Tb_Obs[0, j, i]
    else:
        Tb[k]=Tb_ObsN[0,j,i]
        Error[k]=np.array(Tb_modobs)*Emiss[k]-Tb_ObsN[0,j,i]

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

cmap = plt.get_cmap('rainbow')

'''plt.scatter(Error[:,0],Error[:,1], c=dTdzbar, cmap=cmap)
[plt.text(et-1,em-0.5,site) for et,em,site in zip(Error[:,0], Error[:,1],Sites)]
cbar=plt.colorbar()
cbar.set_label('Mean dT/dz on upper 1000 m', rotation=270, labelpad=15)
plt.xlabel("Error with Tiuri")
plt.ylabel("Error with Mätzler")
plt.grid(which='both')
plt.axis("equal")
plt.xlim(-3,14)
plt.ylim(-3,14)
#plt.autoscale(True)
plt.show()'''

'''plt.bar(dTdzbar, height=Error[:,1], width=5e-4, color='r', label="Mätzler")
plt.bar(dTdzbar-2e-4, height=Error[:,0], width=5e-4, color='b',label="Tiuri")
[plt.text(gradt-1e-3,em+1.5,site, rotation=90) for gradt,em,site in zip(dTdzbar, Error[:,1],Sites)]'''
'''plt.bar(T500_Obs, height=Error[:,1], width=0.75, color='r', label="Mätzler")
plt.bar(T500_Obs-0.25, height=Error[:,0], width=0.75, color='b',label="Tiuri")
[plt.text(gradt-1e-3,em+1.5,site, rotation=90) for gradt,em,site in zip(T500_Obs, Error[:,1],Sites)]'''
plt.grid(which='both')
plt.ylabel("Tb model - Tb SMOS (K), corrected for emissivity")
plt.xlabel("T(d=500 m) (K)")
plt.legend()
#plt.ylim(-2,15)
plt.show()

'''Ts_gr = [ Tz_gr[j,i][0] for i, j in zip(XpixS, YpixS)]
T500_gr = np.array([ np.interp(np.arange(0,H[j,i],100),(np.arange(0,1.05,0.05)*H[j,i]),Tz_gr[j,i])[5] for i, j in zip(XpixS, YpixS)])

Us_Mouginot=np.array([0.1,1,0.1,0.5,3.7,5,0.1])
Us_gr=np.array([0.4,1.4,0.4,8.9,138,16,17])
DeltaUs=Us_gr-Us_Mouginot
DeltaT500_Obs=T500_Obs[:-1]-Ts[:-1]
DeltaT500_gr=T500_gr-Ts_gr
DeltaTs=Ts_gr-Ts[:-1]
DeltaDelta=DeltaT500_gr-DeltaT500_Obs

plt.scatter(DeltaUs, DeltaDelta)
[plt.text(u,t,site) for u,t,site in zip(DeltaUs,DeltaDelta,Sites)]
plt.xlabel("Velocity error (m/a)")
plt.ylabel("Diff")
plt.show()'''

'''plt.scatter(Ts+273.15,Tb, c=Height, cmap=cmap)
[plt.text(ts-1+273.15,tb-3,site) for ts,tb,site in zip(Ts,Tb,Sites)]
cbar=plt.colorbar()
cbar.set_label('Ice thickness (m)', rotation=270, labelpad=15)
plt.xlabel("Ts (K)")
plt.ylabel("Tb SMOS (K)")
plt.plot([200,270],[200,270],'-', c='r', lw=0.5)
plt.grid(which='both')
plt.axis("equal")
plt.xlim(210,260)
plt.ylim(210,260)
#plt.autoscale(True)
plt.show()
#plt.gca().autoscale_view()'''

'''from sklearn import linear_model
regr = linear_model.LinearRegression()
regr.fit(Ts[:,np.newaxis], Tb)
Ts_test = np.linspace(np.min(Ts), np.max(Ts), 100)
plt.plot(Ts_test+273.15, regr.predict(Ts_test[:,np.newaxis]), color='blue', linewidth=1)'''
plt.show()

cmap = plt.get_cmap('inferno')
colors = [cmap(i) for i in Error[:,1]]
#Robin profile
TMcM1=-273.15+np.array([227.742,228.047,228.365,228.696,229.042,229.401,229.776,230.165,230.570,230.990,231.426,231.879,232.347,232.832,233.334,233.853,234.388,234.941,235.512,236.100,236.705,237.328,237.969,238.628,239.304,239.997,240.709,241.437,242.183,242.946,243.726,244.523,245.336,246.165,247.010,247.871,248.747,249.637,250.542,251.461,252.393,253.338,254.296,255.265,256.245,257.236,258.236,259.246,260.265,261.291,262.325,263.365,264.411,265.461,266.516,267.574,268.635,269.698,270.761])
ZMcM1=np.arange(0,50*np.shape(TMcM1)[0],50)
TMcM2=-273.15+np.array([226.469,226.766,227.073,227.392,227.722,228.063,228.417,228.782,229.159,229.549,229.951,230.366,230.794,231.234,231.688,232.155,232.636,233.129,233.637,234.158,234.693,235.241,235.803,236.379,236.969,237.573,238.190,238.822,239.466,240.125,240.797,241.482,242.181,242.893,243.617,244.355,245.105,245.867,246.642,247.428,248.226,249.036,249.856,250.687,251.529,252.380,253.242,254.112,254.991,255.879,256.775,257.678,258.589,259.506,260.430,261.359,262.294,263.233,264.176,265.123,266.074,267.027,267.982,268.939,269.897])
ZMcM2=np.arange(0,50*np.shape(TMcM2)[0],50)

colors=['r','b','g','k','orange','pink','purple',"gray"]

#Plot
i=0
for site, c in zip(Sites, colors):
    #color=colors[i]
    print(site)
    Data = loadtxt("../../SourceData/Temperatures/"+str(site)+".csv", comments="#", delimiter=",",unpack=False)
    Tz_Obs = Data[:, 1]
    depth = Data[:, 0]
    plt.plot(Tz_Obs, depth, label=site, color=c)
    [plt.scatter(Tb_Obs[0,j,i]-273.15, 0,color=c) for i,j,c in zip(XpixS, YpixS, colors)]
    [plt.scatter(Tb_ObsN[0,j,i]-273.15, 0,color="gray") for i,j in zip(XpixN, YpixN)]

    plt.legend()
    i=i+1

'''print((Uxbar[58,153]**2+Uybar[58,153]**2)**0.5)
print((Uxbar[57,149]**2+Uybar[57,149]**2)**0.5)
print((Uxbar[55,146]**2+Uybar[55,146]**2)**0.5)
print((Uxbar[54,142]**2+Uybar[54,142]**2)**0.5)'''

'''print(math.atan(Uybar[58,153]/Uxbar[58,153])/3.14*180)
print(math.atan(Uybar[57,149]/Uxbar[57,149])/3.14*180)
print(math.atan(Uybar[55,146]/Uxbar[55,146])/3.14*180)
print(math.atan(Uybar[54,142]/Uxbar[54,142])/3.14*180)'''
x=[170,171,172,173,174]#77#193#191
y=[49,50,51,52,53]#72#49#55
#x=[160,161,162,163]
#y=[66,66,66,66]
x=[58]#np.arange(202,202,1)
y= [93]#np.arange(90,94,1)
#x=[137,136,135,134]
#y=np.arange(150,145,-1)

#x=[189,189,189,189,189]
#y=[81,82,83,84,85]

'''plt.plot(Tz_gr[58, 153], (1 - Zeta[58,153]) * H[58, 153], c='pink')
plt.plot(Tz_gr[57, 149], (1 - Zeta[57,149]) * H[57, 149], c='r')
plt.plot(Tz_gr[55, 146], (1 - Zeta[55,146]) * H[55, 146], c='g')
plt.plot(Tz_gr[54, 142], (1 - Zeta[54,142]) * H[54, 142], c='b', label='GRISLI')'''
#[plt.plot(Tz_gr[j,i], (1 - Zeta[j,i]) * H[j,i]) for i, j in zip(x,y)]

#[plt.plot(Tz_gr[j,i], (1 - Zeta[j,i]) * H[j,i], '--',color=c) for i, j, k, site, c in zip(XpixS,YpixS, np.arange(0,np.size(XpixS),1),Sites, colors)]
#[print((Uxbar[j,i]**2+Uybar[j,i]**2)**0.5) for i, j in zip(XpixS,YpixS)]

#plt.plot(TMcM1, ZMcM1, c='darkgreen', label="Robin")
#plt.plot(TMcM2, ZMcM2, c='k')
#plt.plot(Tz_gr[78,178],(1-Zeta[78,178])*H[78,178])
#plt.plot(Tz_gr[112,168],(1-Zeta[112,168])*H[112,168])
#plt.plot(Tz_gr[73,79],(1-Zeta[73,79])*H[73,79])
#cbar = plt.colorbar()
#cbar.set_label('Error (K)', rotation=270, labelpad=15)
plt.xlabel("Ice temperature (°C)")
plt.ylabel("Depth (m)")
plt.gca().invert_yaxis()
plt.grid()
plt.legend()
plt.show()

