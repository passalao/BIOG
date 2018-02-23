#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG


Start=time.clock()
#Import GRISLI Data
nc = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
ny=nc.dimensions['y'].size
nx=nc.dimensions['x'].size
nz=nc.dimensions['z'].size

#Extract data from NetCDF file
T = nc.variables['T'][:]  # extract/copy the data
H = nc.variables['H'][:]  # extract/copy the data
T=T[:,:,0:21] #Select layers in the ice, no ground layers wanted.

Tb=np.zeros((nx, ny))
for c in np.arange(0,ny,1):
    if c//BIOG.var.Subsample==float(c)/BIOG.var.Subsample:
        print("y=", c)
        for l in np.arange(0, nx, 1):
            if H[l,c]==1.0:#ocen pixels
                continue
            if l // BIOG.var.Subsample == float(l) / BIOG.var.Subsample:
                print("y=", c, " and x=", l, "time:", time.clock())
                Tz=T[l,c, :]
                Thick=H[l,c]
                Tb[l,c]=BIOG.fun.GetTb_DMRTML(Tz, Thick, BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.NbStreams)

# Export of the enriched GRISLI dataset for KERAS
w_nc_fid = Dataset('../../SourceData/WorkingFiles/GRISLI_Tb_SMOSGrid.nc', 'w', format='NETCDF4')
w_nc_fid.description = "Tb computed from stationary run of GRISLI "
w_nc_fid.createDimension("x", nc.dimensions['x'].size)
w_nc_fid.createDimension("y", nc.dimensions['y'].size)
w_nc_fid.createVariable('Tb','float64',nc.variables['H'].dimensions)
w_nc_fid.variables['Tb'][:] = Tb
w_nc_fid.close()

Stop=time.clock()
print("Elapsed time: ", Stop-Start, 's')