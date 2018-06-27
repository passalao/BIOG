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
from osgeo import gdal

# Import GRISLI data
Model = netCDF4.Dataset('../../SourceData/GRISLI/TB40S004_1_T.nc')
ny_Mod = Model.dimensions['y'].size
nx_Mod = Model.dimensions['x'].size
#nz_Mod = Model.dimensions['z'].size
nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(Model)
Tz=Model.variables['T']
TsRACMO=Tz[0,:,:]

# Import GHF data: already gridded on 40 km
dataset = gdal.Open('../../SourceData/GRISLI/ALBMAP_GhfFM_40km.tif', gdal.GA_ReadOnly)
for x in range(1, dataset.RasterCount + 1):
    band = dataset.GetRasterBand(x)
    GHF = band.ReadAsArray()

#Create new NetCDF file
nc_interp = Dataset('../../SourceData/GRISLI/FoxMaule_40km.nc', 'w', format='NETCDF4')
nc_interp.description = "GHF dataset mapped on GRISLI 40 km grid"
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
nc_interp.createVariable('GHF', 'float64', ('x','y'))
nc_interp.variables['GHF'][:] = GHF[::-1,:]
nc_interp.close()