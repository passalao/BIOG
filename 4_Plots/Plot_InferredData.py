#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

# Import data computed from the model
Model1 = netCDF4.Dataset('../../OutputData/BayesedDataset_Sample10_Bias-5.nc')
ny_Mod1 = Model1.dimensions['y'].size
nx_Mod1 = Model1.dimensions['x'].size
Ts = Model1.variables['Ts']
sigTs = Model1.variables['Ts_std']
PhiG = Model1.variables['PhiG']
sigPhiG = Model1.variables['PhiG_std']
DeltaTb = Model1.variables['DeltaTb']
Cost = Model1.variables['Cost']

#Ts=np.reshape(Ts,(201*225,1))

'''# scatterplot
myplot=plt.scatter(Tb_Obs, Tb_Mod1+offset, c=Acc4Plot, s=1e-3)
cbar=plt.colorbar()
#cbar.set_label('Geothermal flux (mW/m2)', rotation=270, labelpad=15)
cbar.set_label('Accumulation (m/a)', rotation=270, labelpad=15)
plt.plot([0, 270], [0, 270], color='b')
plt.grid()
plt.axis("equal")
plt.autoscale(True)
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb SMOS (K)')
plt.ylabel('Tb GRISLI+'+BIOG.var.RTModel+'(K)')
plt.savefig("../../OutputData/img/Tb_SMOSvsMod_"+BIOG.var.RTModel+".png")
plt.show()'''

'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.Reds
norm = mpl.colors.Normalize(vmin=0, vmax=10)
myplot = plt.pcolormesh(DeltaTb, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0, 10, 1))
cbar.ax.set_xticklabels(['-15', '0', '15'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
#plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()'''

'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.Reds
norm = mpl.colors.Normalize(vmin=0, vmax=20)
myplot = plt.pcolormesh(sigPhiG, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0, 20, 1))
cbar.ax.set_xticklabels(['-15', '0', '15'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
#plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()'''

'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.Reds
norm = mpl.colors.Normalize(vmin=0, vmax=120)
myplot = plt.pcolormesh(PhiG, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0, 120, 10))
cbar.ax.set_xticklabels(['-15', '0', '15'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
#plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()'''