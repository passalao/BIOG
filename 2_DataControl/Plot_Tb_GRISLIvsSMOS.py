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
#Model = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb.nc')
Model = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb_SMOSGrid.nc')
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

'''Lon = np.reshape(Lon, (nx_Obs * ny_Obs, 1))
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
Tb_SMOS = f(xG, yG)  # Interpolate on the new grid'''
offset=0
Tb_Obs=np.array(Tb_Obs)
Tb_Mod=np.array(Tb_Mod)
Error=Tb_Obs[0]-(Tb_Mod+offset)
print(np.shape(Error))

# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.seismic
norm = mpl.colors.Normalize(vmin=-10, vmax=10)
myplot = plt.pcolormesh(Error, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-15, 15, 1))
cbar.ax.set_xticklabels(['-15', '0', '15'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()

'''# scatterplot
Tb_Obs=np.reshape(Tb_Obs,(201*225,1))
Tb_Mod=np.reshape(Tb_Mod,(201*225,1))
plt.scatter(Tb_Obs, Tb_Mod+offset, color='DarkGreen', s=0.001);
plt.plot([0, 270], [0, 270], color='b')
plt.autoscale(True)
plt.text(252,222,"offset = "+str(offset)+" K")
plt.grid()
#plt.axis("equal")
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb SMOS (K)')
plt.ylabel('Tb GRISLI+SMRT (K)')
plt.savefig("../../OutputData/img/Tb_SMOSvsMod_DMRTML.png")
plt.show()'''

plt.close()
