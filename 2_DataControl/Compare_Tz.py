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

GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
H = GRISLI.variables['H']
print(np.shape(H))
Zeta = GRISLI.variables['Zeta']
T = GRISLI.variables['T']

#Dome C coordinates in the SMOS grid
i=161
j=64
#Vostok
i=166
j=83

Tground=T[:,:,20]

print(Tground[150,150])
# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.seismic
norm = mpl.colors.Normalize(vmin=250, vmax=273)
myplot = plt.pcolormesh(Tground+273.15, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(250, 270, 10))
cbar.ax.set_xticklabels(['250', '270'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.show()

Z=Zeta[j,i,:]*H[j,i]
Temp=T[j,i,:]

from numpy import loadtxt
lines = loadtxt("../../SourceData/EDC_Temp_SlantCorrected.csv", comments="#", delimiter=",", unpack=False)

#InterpTemp=np.interp(MonDomaine.z, Temp.Age, Temp.Temp)
#plt.plot(lines[:,1],lines[:,0], label="Measured")
plt.plot(np.array([-57,-53,-48.5,-43,-36])+273.15,3700-np.array([0,500,1000,1500,2000]),label="Measured")
plt.plot(Temp+273.15,Z, label="GRISLI stationnary")
plt.grid()
plt.legend()
plt.xlabel("Temperature (K)")
plt.ylabel("Z (m)")
plt.savefig("../../OutputData/img/Tz_at_Vostok.png")
plt.show()