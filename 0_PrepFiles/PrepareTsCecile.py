#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import NC_Resources as ncr
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.interpolate as si

TsMAR = netCDF4.Dataset('../../SourceData/MAR/MAR-ERA-Interim_ST_20c_ave.nc4')
Ts=TsMAR.variables['ST']
nx = TsMAR.dimensions['y'].size
ny = TsMAR.dimensions['x'].size

XMAR=np.arange(-(nx-1)/2*35000, ((nx-1)/2+1)*35000, 35000)
YMAR=np.arange(-(ny-1)/2*35000, ((ny-1)/2+1)*35000, 35000)

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']

#Interpolation on SMOS Grid
f = si.interp2d(YMAR, XMAR, Ts, kind='cubic')
TsInterp = f(X[0, :], Y[:, 0])
TsInterp=TsInterp*(4-np.array(Mask))/3

#Write output NetCDF file
outfile = r'../../SourceData/WorkingFiles/TsMAR.nc'
cols = len(X[0,:])
rows = len(Y[:,0])
dsout = netCDF4.Dataset(outfile, 'w', clobber=True)
Yout = dsout.createDimension('y', rows)
Yout = dsout.createVariable('y', 'f4', ('y',))
Yout.standard_name = 'y'
Yout.units = 'm'
Yout.axis = "Y"
Yout[:] = Y[:,0]+25000
Xout = dsout.createDimension('x', cols)
Xout = dsout.createVariable('x', 'f4', ('x',))
Xout.standard_name = 'x'
Xout.units = 'm'
Xout.axis = "X"
Xout[:] = X[0,:]-25000

dsout.createVariable('TsMAR','float64',('y','x'))
dsout.variables['TsMAR'][:] = np.array(TsInterp[::-1,:])
crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'''
