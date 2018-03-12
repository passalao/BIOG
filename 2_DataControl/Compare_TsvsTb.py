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
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
Lon = Obs.variables['lon']
Lat = Obs.variables['lat']
Mask = Obs.variables['mask']

Obs2 = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
Tb_Obs2 = Obs2.variables['BT_V']


#Import GRISLI Data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
Acc = GRISLI.variables['Acc']
T = GRISLI.variables['T']
Ts=T[:,:,0]

#Import Tb data out of GRISLI
Model1 = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb_SMOSGrid_'+BIOG.var.RTModel+'_'+'Matzler'+'.nc')
ny_Mod1 = Model1.dimensions['y'].size
nx_Mod1 = Model1.dimensions['x'].size
Tb_Mod1 = Model1.variables['Tb']

#Select data in the continent
Tb_Obs=Tb_Obs*(4-np.array(Mask))/3
Tb_Obs2=Tb_Obs2*(4-np.array(Mask))/3
Ts=(Ts+273.15)*(4-np.array(Mask))/3
Tb_Mod1=Tb_Mod1*(4-np.array(Mask))/3
DeltaT=Tb_Obs-Ts
#DeltaT=Tb_Obs-Tb_Obs2

Tb_Obs=np.array(Tb_Obs)
Tb_Mod1=np.array(Tb_Mod1)
Ts=np.array(Ts)

Tb_Mod1=np.reshape(Tb_Mod1,(201*225,1))
Tb_Obs=np.reshape(Tb_Obs,(201*225,1))
Ts=np.reshape(Ts,(201*225,1))

# scatterplot
myplot=plt.scatter(Tb_Mod1, Ts, c="Red", s=0.01)
#cbar=plt.colorbar()
#cbar.set_label('Geothermal flux (mW/m2)', rotation=270, labelpad=15)
#cbar.set_label('Accumulation (m/a)', rotation=270, labelpad=15)
plt.plot([0, 270], [0, 270], color='b')
plt.grid()
plt.axis("equal")
plt.autoscale(True)
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb DMRT-ML (K)')
plt.ylabel('Ts GRISLI (K)')
#plt.savefig("../../OutputData/img/TbSMRTTiuri_vs_TsGRISLI.png")
plt.show()

# scatterplot
myplot=plt.scatter(Tb_Mod1, Ts, s=1e-1)
#cbar=plt.colorbar()
#cbar.set_label('Geothermal flux (mW/m2)', rotation=270, labelpad=15)
#cbar.set_label('Accumulation (m/a)', rotation=270, labelpad=15)
plt.plot([0, 270], [0, 270], color='b')
plt.grid()
plt.axis("equal")
plt.autoscale(True)
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb SMOS (K)')
plt.ylabel('Ts GRISLI (K)')
#plt.savefig("../../OutputData/img/TbSMOS_vs_TsGRISLI.png")
plt.show()

'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.seismic#Reds_r#seismic
norm = mpl.colors.Normalize(vmin=-20, vmax=20)
myplot = plt.pcolormesh(DeltaT[0,:,:], cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-40, 1, 5))
cbar.ax.set_xticklabels(['-15', '0', '15'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
#plt.savefig("../../OutputData/img/Delta_TbSMOS-TsGRISLI.png")
plt.show()'''

plt.close()