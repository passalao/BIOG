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
import pyproj
import scipy.interpolate as si
import NC_Resources as ncr

# Import GRISLI data
Model = netCDF4.Dataset('../../SourceData/GRISLI/TB40S004_1_T.nc')
ny_Mod = Model.dimensions['y'].size
nx_Mod = Model.dimensions['x'].size
#nz_Mod = Model.dimensions['z'].size
nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(Model)
Tz=Model.variables['T']
TsRACMO=Tz[0,:,:]

# Import Crocus data
ObsCrocus = netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc')
Ts = ObsCrocus.variables['TsCrocus'][::-1,:]
TsCrocus=Ts-273.15
TsCrocus[TsCrocus==-273.15]=0

# Import SMOS data for grid coordinates
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
Lon = Obs.variables['lon']
Lat = Obs.variables['lat']
Mask = Obs.variables['mask']
x=Obs.variables['x_ease2']
y=Obs.variables['y_ease2']

nc_obsattrs, nc_obsdims, nc_obsvars = ncr.ncdump(Obs)

Lon = np.reshape(Lon, (nx_Obs * ny_Obs, 1))
Lat = np.reshape(Lat, (nx_Obs * ny_Obs, 1))

# GRISLI coordinates
DeltaX=0.5*(ny_Mod*40e3-ny_Obs*15e3)
DeltaY=0.5*(nx_Mod*40e3-nx_Obs*15e3)
print(DeltaX, DeltaY)
xG = 40e3*np.linspace(-(nx_Mod-1)/2,(nx_Mod-1)/2, nx_Mod)
yG = 40e3*np.linspace((ny_Mod-1)/2,-(ny_Mod-1)/2 , ny_Mod)

# SMOS coordinates
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj("+init=EPSG:3031")  #
XX, YY = pyproj.transform(wgs84, StereoPol, Lon, Lat)
X = np.reshape(XX, (nx_Obs, ny_Obs))
Y = np.reshape(YY, (nx_Obs, ny_Obs))
Tb_Obs=Tb_Obs[0,:,:]

#Create new NetCDF file for interpolated data
nc_interp = Dataset('../../SourceData/WorkingFiles/TsMappedon40km.nc', 'w', format='NETCDF4')
nc_interp.description = "GRISLI data mapped on SMOS grid"
print(Model.dimensions["x"].size)
#Create NetCDF dimensions
for dim in nc_moddims:
    dim_name=dim
    if dim=="rows":
        dim_name="x"
    if dim=="cols":
        dim_name="y"
    nc_interp.createDimension(dim_name, Model.dimensions[dim_name].size)

#Interpolates Ts on GRISLI grid
nc_interp.createVariable('TsCrocus', 'float64', ('x','y'))
f = si.interp2d(X[0, :], Y[:, 0], TsCrocus, kind='cubic')
Interp_TsCrocus = f(xG, yG)[::-1, :]

Mask=np.ones(np.shape(Interp_TsCrocus))
Mask[Interp_TsCrocus>-8]=0
Interp_TsCrocus=Mask*Interp_TsCrocus+(1-Mask)*TsRACMO

nc_interp.variables['TsCrocus'][:] = Interp_TsCrocus

# Geographic plot
fig, ax = plt.subplots(nrows=1, ncols=1, sharex='col')
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=-70, vmax=0)
myplot=ax.pcolormesh(Interp_TsCrocus, cmap=cmap, norm=norm)

plt.show()

'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=180, vmax=250)
myplot = plt.pcolormesh(Tb_Obs, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-10, 10, 1))
cbar.ax.set_xticklabels(['-10', '0', '10'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()
plt.close()'''
nc_interp.close()
