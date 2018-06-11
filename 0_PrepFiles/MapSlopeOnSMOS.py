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

#Import Slopes
Slopes=netCDF4.Dataset('../../SourceData/WorkingFiles/TsCrocusModulations_and_H.nc')
ny_Mod = Slopes.dimensions['y'].size
nx_Mod = Slopes.dimensions['x'].size
GradS=Slopes.variables['Slope']
nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(Slopes)

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
DeltaX=0.5*(ny_Mod*25e3-ny_Obs*15e3)
DeltaY=0.5*(nx_Mod*25e3-nx_Obs*15e3)
print(DeltaX, DeltaY)
xG = 25e3*np.linspace(-(nx_Mod-1)/2,(nx_Mod-1)/2, nx_Mod)
yG = 25e3*np.linspace((ny_Mod-1)/2,-(ny_Mod-1)/2 , ny_Mod)

# SMOS coordinates
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj("+init=EPSG:3031")  #
XX, YY = pyproj.transform(wgs84, StereoPol, Lon, Lat)
X = np.reshape(XX, (nx_Obs, ny_Obs))
Y = np.reshape(YY, (nx_Obs, ny_Obs))

#Create new NetCDF file for interpolated data
nc_interp = Dataset('../../SourceData/WorkingFiles/SurfaceSlopesOnSMOS.nc', 'w', format='NETCDF4')
nc_interp.description = "Slopes data mapped on SMOS grid"
print(Slopes.dimensions["x"].size)

#Create NetCDF dimensions
for dim in nc_obsdims:
    dim_name=dim
    if dim=="rows":
        dim_name="x"
    if dim=="cols":
        dim_name="y"
    nc_interp.createDimension(dim_name, Obs.dimensions[dim].size)

#Interpolates Slope on SMOS grid
nc_interp.createVariable('Slopes', 'float64', ('x','y'))
f = si.interp2d(xG,yG, GradS, kind='cubic')
Interp_GradS = f(X[0, :], Y[:, 0])

nc_interp.variables['Slopes'][:] = Interp_GradS

# Geographic plot
fig, ax = plt.subplots(nrows=1, ncols=1, sharex='col')
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=-1e-0, vmax=1e-0)
myplot=ax.pcolormesh(Interp_GradS, cmap=cmap, norm=norm)
plt.show()

nc_interp.close()
