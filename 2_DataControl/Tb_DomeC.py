#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
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
Mask = Obs.variables['mask']

#Import GRISLI data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
H = GRISLI.variables['H']
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']

f = open('RT_Sites.csv', 'w')

f.write('Number of vertical layers: '+str(BIOG.var.NbLayers)+"\n")
f.write('Permittivity model: '+str(BIOG.var.Perm)+"\n")
f.write('Number of streams: '+str(BIOG.var.Angle)+"\n")
f.write('Site Tb Ts_LocalObs Ts_Dataset Tb_LocalObs Tb_Dataset \n')

#Dome C
print("         Dome C         ")
i=65
j=167
print("Tb SMOS at Dome C:", Tb_Obs[0,i,j])
Tz_gr_at_Point=Tz_gr[i,j,:]
print("Ts:", Tz_gr_at_Point[0]+273.15)
Data = loadtxt("../../SourceData/Temperatures/EDC_Temp_SlantCorrected.csv", comments="#", delimiter=",", unpack=False)
Tz_Obs = Data[:,1]-273.15
Z=Data[:,0]
Tb_mod = BIOG.fun.GetTb(Tz_Obs[::-1], Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz GRISLI: ", Tb_gr)
print("Biais d� � GRISLI :", Tb_gr-Tb_mod)
print("Biais d� � DMRT-ML :", Tb_mod-Tb_Obs[0,i,j])
print(' ')
f.write("Dome C "+str(Tb_Obs[0,i,j])+" "+str(Tz_Obs[-1]+273.15)+" "+str(Tz_gr_at_Point[0]+273.15)+" "+str(Tb_Obs[0,i,j])+" "+str(Tb_gr)+"\n")

#Vostok
print("         Vostok         ")
i=86
j=165
print("Tb SMOS at Vostok:", Tb_Obs[0,i,j])
Tz_gr_at_Point=Tz_gr[i,j,:]
Tz_Obs=[-57,-53,-48.5,-43,-36,-2]
Z=3450.-np.array([0,500,1000,1500,2000,3450])
Tb_mod = BIOG.fun.GetTb(Tz_Obs, Z[0],BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[0], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz GRISLI: ", Tb_gr)
print("Biais d� � GRISLI :", Tb_gr-Tb_mod)
print("Biais d� � DMRT-ML :", Tb_mod-Tb_Obs[0,i,j])
print(' ')
f.write("Vostok "+str(Tb_Obs[0,i,j])+" "+str(Tz_Obs[0]+273.15)+" "+str(Tz_gr_at_Point[0]+273.15)+" "+str(Tb_Obs[0,i,j])+" "+str(Tb_gr)+"\n")

#Dome Fuji
print("         Dome Fuji         ")
i=143
j=147
print("Tb SMOS at Dome Fuji:", Tb_Obs[0,i,j])
Tz_gr_at_Point=Tz_gr[i,j,:]
Data = loadtxt("../../SourceData/Temperatures/DF_Temp.csv", comments="#", delimiter=" ", unpack=False)
Tz_Obs=Data[:,1]
Z=Data[:,0]
#plt.plot(Tz_Obs,Z[-1]-Z, color='r')
#plt.plot(Tz_gr_at_Point, Zeta[143,147,:]*Z[-1], color='b')
#plt.show()
Tb_mod = BIOG.fun.GetTb(Tz_Obs, Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz GRISLI: ", Tb_gr)
print("Biais d� � GRISLI :", Tb_gr-Tb_mod)
print("Biais d� � DMRT-ML :", Tb_mod-Tb_Obs[0,i,j])
print(' ')
f.write("Dome Fuji "+str(Tb_Obs[0,i,j])+" "+str(Tz_Obs[0]+273.15)+" "+str(Tz_gr_at_Point[0]+273.15)+" "+str(Tb_Obs[0,i,j])+" "+str(Tb_gr)+"\n")


print("         Law Dome         ")
i=62
j=206
print("Tb SMOS at Law Dome:", Tb_Obs[0,i,j])
Tz_gr_at_Point=Tz_gr[i,j,:]
print("Ts:", Tz_gr_at_Point[0]+273.15)
Data = loadtxt("../../SourceData/Temperatures/LawDome_Temp.csv", comments="#", delimiter=" ", unpack=False)
Tz_Obs=Data[:,1]
Z=Data[:,0]
Z=Z[::-1]
#plt.plot(Tz_Obs,Z, color='r')
#plt.plot(Tz_gr_at_Point, Zeta[i,j,:]*Z[-1], color='b')
#plt.show()
Tb_mod = BIOG.fun.GetTb(Tz_Obs[::-1], Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz GRISLI: ", Tb_gr)
print("Biais d� � GRISLI :", Tb_gr-Tb_mod)
print("Biais d� � DMRT-ML :", Tb_mod-Tb_Obs[0,i,j])
print(' ')

print("         T-UC         ")
i=80
j=94
print("Tb SMOS at T-UC:", Tb_Obs[0,i,j])
Tz_gr_at_Point=Tz_gr[i,j,:]
print("Ts:", Tz_gr_at_Point[0]+273.15)
Data = loadtxt("../../SourceData/Temperatures/T-UC-1993-14.txt", comments="#", delimiter=" ", unpack=False)
Tz_Obs=Data[:,0]
Z=Data[:,1]
#plt.plot(Tz_Obs,Z, color='r')
#plt.plot(Tz_gr_at_Point, Zeta[80,94,:]*Z[-1], color='b')
#plt.show()
Tb_mod = BIOG.fun.GetTb(Tz_Obs[::-1], Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz GRISLI: ", Tb_gr)
print("Biais d� � GRISLI :", Tb_gr-Tb_mod)
print("Biais d� � DMRT-ML :", Tb_mod-Tb_Obs[0,i,j])
print(' ')

print("         T-AIS         ")
i=79
j=93
print("Tb SMOS at T-AIS:", Tb_Obs[0,i,j])
Tz_gr_at_Point=Tz_gr[i,j,:]
print("Ts:", Tz_gr_at_Point[0]+273.15)
Data = loadtxt("../../SourceData/Temperatures/T-UC-1993-14.txt", comments="#", delimiter=" ", unpack=False)
Tz_Obs=Data[:,0]
Z=Data[:,1]
#plt.plot(Tz_Obs,Z, color='r')
#plt.plot(Tz_gr_at_Point, Zeta[79,93,:]*Z[-1], color='b')
#plt.show()
Tb_mod = BIOG.fun.GetTb(Tz_Obs[::-1], Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
print("Tb mod�le avec Tz GRISLI: ", Tb_gr)
print("Biais d� � GRISLI :", Tb_gr-Tb_mod)
print("Biais d� � DMRT-ML :", Tb_mod-Tb_Obs[0,i,j])
print(' ')

'''#Fictitious point, for verification
print("         Fictitious point         ")
T=220-273.15
H=3000
print("Tb SMOS :", T+273.15)
Tz_gr_at_Point=T*np.ones(10)
angles=np.linspace(50,61,21)
Tb_gr=np.zeros(21)
for a in angles:
    Tb_gr[np.where(angles==a)[0]] = BIOG.fun.GetTb(Tz_gr_at_Point, H, 10, 1.4e9, a, 64, "Tiuri", "DMRT-ML")
plt.plot(angles, Tb_gr)
plt.show()
#Tb_gr= BIOG.fun.GetTb(Tz_gr_at_Point, H, BIOG.var.NbLayers, BIOG.var.Freq, 60.6, 256, BIOG.var.Perm, BIOG.var.RTModel)
#print("Tb mod�le avec Tz GRISLI (beware of density): ", Tb_gr)
print(' ')'''


'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.seismic
norm = mpl.colors.Normalize(vmin=250, vmax=273)
myplot = plt.pcolormesh(Tground+273.15, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(250, 270, 10))
cbar.ax.set_xticklabels(['250', '270'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.show()'''

#InterpTemp=np.interp(MonDomaine.z, Temp.Age, Temp.Temp)
'''#plt.plot(lines[:,1],lines[:,0], label="Measured")
plt.plot(np.array([-57,-53,-48.5,-43,-36])+273.15,3700-np.array([0,500,1000,1500,2000]),label="Measured")
plt.plot(Temp+273.15,Z, label="GRISLI stationnary")
plt.grid()
plt.legend()
plt.xlabel("Temperature (K)")
plt.ylabel("Z (m)")
#plt.savefig("../../OutputData/img/Tz_at_Vostok.png")
plt.show()'''
f.close()