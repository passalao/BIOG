#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
from netCDF4 import Dataset
import netCDF4, BIOG
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
import scipy.interpolate as si
import NC_Resources as ncr

# Import GRISLI data
#Model = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S004_1_KERAS.nc')
Model = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/Corrected_Tz_GRISLI.nc')
ny_Mod = Model.dimensions['y'].size
nx_Mod = Model.dimensions['x'].size
nz_Mod = Model.dimensions['z'].size
nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(Model)

# Import SMOS data
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

'''xG = np.linspace(-2.805e6, 3e6, nx_Mod)
yG = np.linspace(2.805e6, -2.805e6, ny_Mod)'''

# SMOS coordinates
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj("+init=EPSG:3031")  #
XX, YY = pyproj.transform(wgs84, StereoPol, Lon, Lat)
X = np.reshape(XX, (nx_Obs, ny_Obs))
Y = np.reshape(YY, (nx_Obs, ny_Obs))
Tb_Obs=Tb_Obs[0,:,:]

#Create new NetCDF fil for interpolated data
nc_interp = Dataset('../../SourceData/GRISLI/Avec_FoxMaule/Corrected_Tz_MappedonSMOS.nc', 'w', format='NETCDF4')
nc_interp.description = "GRISLI data mapped on SMOS grid"

#Create NetCDF dimensions
for dim in nc_obsdims:
    dim_name=dim
    if dim=="rows":
        dim_name="x"
    if dim=="cols":
        dim_name="y"
    nc_interp.createDimension(dim_name, Obs.dimensions[dim].size)
nc_interp.createDimension('z', Model.dimensions['z'].size)

#Interpolates GRISLI on SMOS grid
for var in Model.variables:
    print(var)
    if Model.variables[var].dimensions==('y','x'):
        print("Interpolating ", var)
        nc_interp.createVariable(var, 'float64', ('x','y'))
        print(np.shape(xG),np.shape(xG),np.shape(Model.variables[var][:]))
        f = si.interp2d(xG, yG, Model.variables[var][:], kind='cubic')
        Interp_Var = f(X[0,:],Y[:,0])
        nc_interp.variables[var][:] = Interp_Var[::-1,:]
    if Model.variables[var].dimensions==('z','y','x') or Model.variables[var].dimensions==('zm','y','x'):
        print("Interpolating ", var)
        for k in np.arange(nz_Mod):
            if k==0:
                nc_interp.createVariable(var, 'float64', ('x','y','z'))
            f = si.interp2d(xG, yG, Model.variables[var][k,:,:], kind='cubic')
            Interp_Var = f(X[0,:],Y[:,0])
            nc_interp.variables[var][:,:,k] = Interp_Var[::-1,:]
nc_interp.close()

'''# Geographic plot
fig, ax = plt.subplots(nrows=2, ncols=2, sharex='col')
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=180, vmax=250)
myplot=ax[0,0].pcolormesh(Tb_Obs, cmap=cmap, norm=norm)
#cbar = fig.colorbar(myplot, ticks=np.arange(-10, 10, 1))
#cbar.ax.set_xticklabels(['-10', '0', '10'])  # vertically oriented colorbar
norm = mpl.colors.Normalize(vmin=0, vmax=3500)
ax[1,0].pcolormesh(nc_interp.variables['S'], cmap=cmap, norm=norm)
ax[1,1].pcolormesh(Model.variables['S'], cmap=cmap, norm=norm)
plt.autoscale(True)
plt.axis('equal')
ax[0, 0].set_xlim([0, 225])
ax[0, 0].set_ylim([0, 201.])
ax[1, 1].set_xlim([0, 141.])
ax[1, 1].set_ylim([0, 201.])
#plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()'''

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
