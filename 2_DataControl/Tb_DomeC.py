#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
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

#Dome C
print("         Dome C         ")
print("Tb SMOS at Dome C:", Tb_Obs[0,65,167])
Tz_gr_at_Point=Tz_gr[65,167,:]
Data = loadtxt("../../SourceData/EDC_Temp_SlantCorrected.csv", comments="#", delimiter=",", unpack=False)
Tz_Obs=Data[:,1]-273.15
Z=Data[:,0]
Tb_mod = BIOG.fun.GetTb(Tz_Obs[::-1], Z[-1], 10, 1.4e9, 52.5, 16, "Tiuri", "DMRT-ML")
print("Tb modèle avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], 10, 1.4e9, 52.5, 16, "Tiuri", "DMRT-ML")
print("Tb modèle avec Tz GRISLI: ", Tb_gr)
print("Biais du à GRISLI :", Tb_gr-Tb_mod)
print(' ')

#Vostok
print("         Vostok         ")
print("Tb SMOS at Vostok:", Tb_Obs[0,86,165])
Tz_gr_at_Point=Tz_gr[86,165,:]
Tz_Obs=[-57,-53,-48.5,-43,-36,-2]
Z=3450.-np.array([0,500,1000,1500,2000,3450])
Tb_mod = BIOG.fun.GetTb(Tz_Obs, Z[-1], 10, 1.4e9, 52.5, 16, "Tiuri", "DMRT-ML")
print("Tb modèle avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], 10, 1.4e9, 52.5, 16, "Tiuri", "DMRT-ML")
print("Tb modèle avec Tz GRISLI: ", Tb_gr)
print("Biais du à GRISLI :", Tb_gr-Tb_mod)
print(' ')

#Dome Fuji
print("         Dome Fuji         ")
print("Tb SMOS at Dome Fuji:", Tb_Obs[0,143,147])
Tz_gr_at_Point=Tz_gr[143,147,:]
Data = loadtxt("../../SourceData/DF_Temp.csv", comments="#", delimiter=" ", unpack=False)
Tz_Obs=Data[:,1]
Z=Data[:,0]
plt.plot(Tz_Obs,Z[-1]-Z, color='r')
plt.plot(Tz_gr_at_Point, Zeta[143,147,:]*Z[-1], color='b')
plt.show()
Tb_mod = BIOG.fun.GetTb(Tz_Obs, Z[-1], 10, 1.4e9, 52.5, 16, "Tiuri", "DMRT-ML")
print("Tb modèle avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], 10, 1.4e9, 52.5, 16, "Tiuri", "DMRT-ML")
print("Tb modèle avec Tz GRISLI: ", Tb_gr)
print("Biais du à GRISLI :", Tb_gr-Tb_mod)
print(' ')

#Byrd
print("         Byrd         ")
print("Tb SMOS at Byrd:", Tb_Obs[0,80,74])
Tz_gr_at_Point=Tz_gr[80,74,:]
Data = [[0,-28.3],[690,-22.]]
Tz_Obs=[-28.3,-22]
Z=[0,690]
Tb_mod = BIOG.fun.GetTb(Tz_Obs, Z[-1], 10, 1.4e9, 52.5, 16, "Tiuri", "DMRT-ML")
print("Tb modèle avec Tz obs: ", Tb_mod)
Tb_gr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], 10, 1.4e9, 52.5, 16, "Tiuri", "DMRT-ML")
print("Tb modèle avec Tz GRISLI: ", Tb_gr)
print("Biais du à GRISLI :", Tb_gr-Tb_mod)
print(' ')



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