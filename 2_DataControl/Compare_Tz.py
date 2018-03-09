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
Zeta = GRISLI.variables['Zeta']
T = GRISLI.variables['T']
#Z=np.array([np.array(H).T]*21).T*np.array(Zeta)

rho=917
g=9.81
PMP=273.16-0.074*rho*g*(H[:,:]-30)/1e6-0.024
print(PMP[150,150])
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
cmap = mpl.cm.seismic#Reds_r
norm = mpl.colors.Normalize(vmin=-10, vmax=10)
myplot = plt.pcolormesh(Tground+273.15-PMP, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-10, 10, 5))
cbar.ax.set_xticklabels(['250', '270'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.show()

Z=Zeta[j,i,:]*H[j,i]
Temp=T[j,i,:]

from numpy import loadtxt
lines = loadtxt("../../SourceData/EDC_Temp_SlantCorrected.csv", comments="#", delimiter=",", unpack=False)

'''#InterpTemp=np.interp(MonDomaine.z, Temp.Age, Temp.Temp)
#plt.plot(lines[:,1],lines[:,0], label="Measured")
plt.plot(np.array([-57,-53,-48.5,-43,-36])+273.15,3700-np.array([0,500,1000,1500,2000]),label="Measured")
plt.plot(Temp+273.15,Z, label="GRISLI stationnary")
plt.grid()
plt.legend()
plt.xlabel("Temperature (K)")
plt.ylabel("Z (m)")
#plt.savefig("../../OutputData/img/Tz_at_Vostok.png")
plt.show()'''