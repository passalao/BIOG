#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy.interpolate as si
import NC_Resources as ncr
import seaborn as sns
import pandas as pd

# Import Tb data computed from the model
#Model = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb.nc')
Model = netCDF4.Dataset('../../SourceData/GRISLI/RGL-0194_an2040.nc')
ny = Model.dimensions['y'].size
nx = Model.dimensions['x'].size
Ts_out = Model.variables['T']

SourceData = netCDF4.Dataset('../../SourceData/GRISLI/temp_15km.grd')
Ts_in=SourceData.variables['z']

Ecart=np.array(Ts_out[0,:,:])-np.array(Ts_in)

# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.seismic
norm = mpl.colors.Normalize(vmin=-5, vmax=5)
myplot = plt.pcolormesh(Ecart, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-5, 5, 1))
cbar.ax.set_xticklabels(['-5', '0', '5'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.savefig("../../OutputData/img/Ecart_TsinTsout.png")
#plt.show()

# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.seismic
norm = mpl.colors.Normalize(vmin=-55, vmax=0)
myplot = plt.pcolormesh(Ts_in, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-55, 0, 5))
cbar.ax.set_xticklabels(['-55', '0'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.show()