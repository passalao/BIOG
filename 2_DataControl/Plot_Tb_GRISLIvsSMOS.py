#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy.interpolate as si
import NC_Resources as ncr

# Import GRISLI data
Model = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb.nc')
ny_Mod = Model.dimensions['y'].size
nx_Mod = Model.dimensions['x'].size
Tb_Mod = Model.variables['Tb']

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
f = si.interp2d(X[0, :], Y[:, 0], Tb_Obs, kind='cubic')
Tb_SMOS = f(xG, yG)  # Interpolate on the new grid

# scatterplot
import seaborn
#x=np.reshape(Tb_SMOS.T[::-1], (nx_Mod*ny_Mod, 1))
#y=np.reshape(Tb_Mod, (nx_Mod*ny_Mod, 1))
#seaborn.kdeplot(x[:,0],y[:,0])
plt.scatter(Tb_SMOS[::-1], Tb_Mod, color='DarkGreen', s=1);
plt.plot([0, 270], [0, 270], color='b')
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb SMOS (K)')
plt.ylabel('Tb GRISLI+SMRT (K)')
plt.savefig("../../OutputData/img/Tb_SMOSvsMod.png")
plt.show()

'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=200, vmax=270)
myplot = plt.pcolormesh(Tb_SMOS[::-1], cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(200, 270, 10))
cbar.ax.set_xticklabels(['200', '250'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.show()'''
