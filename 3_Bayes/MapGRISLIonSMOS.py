#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, BIOG
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy.interpolate as si
import NC_Resources as ncr

# Import GRISLI data
Model = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI4KERAS.nc')
ny_Mod = Model.dimensions['y'].size
nx_Mod = Model.dimensions['x'].size
#for var in Model.variables:
PhiG= Model.variables['PhiG']

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
Lon = Obs.variables['lon']
Lat = Obs.variables['lat']
Mask = Obs.variables['mask']

Lon = np.reshape(Lon, (nx_Obs * ny_Obs, 1))
Lat = np.reshape(Lat, (nx_Obs * ny_Obs, 1))

# GRISLI coordinates
xG = np.linspace(-2.805e6, 3e6, 387)
yG = np.linspace(2.805e6, -2.805e6, 374)

# SMOS coordinates
import mpl_toolkits.basemap.pyproj as pyproj
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj("+init=EPSG:3031")  #
XX, YY = pyproj.transform(wgs84, StereoPol, Lon, Lat)
X = np.reshape(XX, (nx_Obs, ny_Obs))
Y = np.reshape(YY, (nx_Obs, ny_Obs))
Tb_Obs=Tb_Obs[0,:,:]

#Pour interpoler GRISLI sur SMOS
f = si.interp2d(xG, yG, PhiG, kind='cubic')
PhiG_interp = f(Y[0,:],X[:,0])  # Interpolate on the SMOS grid
print(np.shape(PhiG_interp), np.shape(Tb_Obs), np.shape(X), np.shape(Y))

# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=180, vmax=250)
myplot = plt.pcolormesh(Tb_Obs, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-10, 10, 1))
cbar.ax.set_xticklabels(['-10', '0', '10'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
#plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()

norm = mpl.colors.Normalize(vmin=40, vmax=120)
myplot = plt.pcolormesh(PhiG, cmap=cmap, norm=norm)
plt.autoscale(True)
plt.axis('equal')
plt.show()
plt.close()