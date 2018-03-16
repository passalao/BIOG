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

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
Mask = Obs.variables['mask']

Start=time.time()
#Import GRISLI Data
nc = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
ny=nc.dimensions['y'].size
nx=nc.dimensions['x'].size
nz=nc.dimensions['z'].size

#Extract data from NetCDF file
T = nc.variables['T'][:]  # extract/copy the data
H = nc.variables['H'][:]  # extract/copy the data
T=T[:,:,0:21] #Select layers in the ice, no ground layers wanted.

Tb1=np.zeros((nx, ny))
Tb2=np.zeros((nx, ny))

for c in np.arange(0,ny,1):
    if c//BIOG.var.Subsample==float(c)/BIOG.var.Subsample:
        for l in np.arange(0, nx, 1):
            if Mask[l,c]==1.0 and l // BIOG.var.Subsample == float(l) / BIOG.var.Subsample:
                print("y=", c, " and x=", l, "time:", time.time()-Start)
                Tz=T[l,c, :]
                Thick=H[l,c]
                Tb1[l,c]=BIOG.fun.GetTb(Tz, Thick, BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
                Tb2[l,c]=BIOG.fun.GetTb(Tz-1, Thick, BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams, BIOG.var.Perm, BIOG.var.RTModel)
                #Tb[l,c]=BIOG.fun.GetTb_SMRT(Tz, Thick, BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.Perm)
                if Tb1[l,c]-Tb2[l,c]<0.99:
                    print("Emissivity not 1", )

'''plt.scatter(np.reshape(Tb,(201*225,1)), np.reshape(Tb2,(201*225,1)),s=1)
plt.plot([0,300],[0,300], c="r")
plt.xlim(200,270)
plt.ylim(200,270)
plt.show()'''

# Export of the enriched GRISLI dataset for KERAS
#w_nc_fid = Dataset('../../SourceData/WorkingFiles/GRISLI_Tb_SMOSGrid_'+BIOG.var.RTModel+'_'+BIOG.var.Perm+'_test.nc', 'w', format='NETCDF4')
w_nc_fid = Dataset('../../SourceData/WorkingFiles/TbMod_and_Emissivity_test.nc', 'w', format='NETCDF4')
w_nc_fid.description = "Tb computed from stationary run of GRISLI with "+ BIOG.var.RTModel
w_nc_fid.createDimension("x", nc.dimensions['x'].size)
w_nc_fid.createDimension("y", nc.dimensions['y'].size)
w_nc_fid.createVariable('Tb','float64',nc.variables['H'].dimensions)
w_nc_fid.variables['Tb'][:] = Tb1
w_nc_fid.createVariable('Emissivity','float64',nc.variables['H'].dimensions)
w_nc_fid.variables['Emissivity'][:] = Tb2-Tb1
w_nc_fid.close()

Stop=time.time()
print("Elapsed time: ", Stop-Start, 's')