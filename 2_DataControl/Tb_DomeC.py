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

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
Lon = Obs.variables['lon']
Lat = Obs.variables['lat']
Mask = Obs.variables['mask']

print("Tb at Dome C:", Tb_Obs[0,64,161])

#Import GRISLI data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
H = GRISLI.variables['H']
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Tz_gr_DomeC=Tz_gr[64,161,:]

from numpy import loadtxt
DomeC_Data = loadtxt("../../SourceData/EDC_Temp_SlantCorrected.csv", comments="#", delimiter=",", unpack=False)
Tz_Obs=DomeC_Data[:,1]-273.15
Z=DomeC_Data[:,0]

Tb_mod = BIOG.fun.GetTb(Tz_Obs[::-1], Z[-1], 10, 1.4e9, 52.5, 16, "Matzler", "DMRT-ML")
print("Tb modèle avec Tz obs: ", Tb_mod)

Tb_gr = BIOG.fun.GetTb(Tz_gr_DomeC, Z[-1], 10, 1.4e9, 52.5, 16, "Matzler", "DMRT-ML")
print("Tb modèle avec Tz GRISLI: ", Tb_gr)

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