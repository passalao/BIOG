from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
from io import StringIO

Transect="DC2Triangle"#DC2Triangle" #DC2McM, DC2VL or DC2Triangle

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
Tb_Obs = Obs.variables['BT_V']
X=Obs.variables['x_ease2']
Y=Obs.variables['y_ease2']
Lx=25e3
Ly=25e3

#Import transect coordinates
TransCoord=np.loadtxt("../../SourceData/Macelloni2016/Macelloni2016_DataTransects/"+Transect+".csv", comments="#", delimiter=" ")
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj(init="EPSG:6932")
Xs, Ys = pyproj.transform(wgs84, StereoPol, TransCoord[:,0], TransCoord[:,1])
Xpix=(Xs-X[0][0])//Lx+1
Ypix=(Ys-Y[-1][-1])//Ly+1
Ts=np.zeros(np.shape(TransCoord)[0])

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles_GRISLIini/GRISLIMappedonSMOS.nc')
H = GRISLI.variables['H']
S = GRISLI.variables['S']
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
H_trans=[H[j, i] for i,j in zip(Xpix,Ypix)]
S_trans=[S[j, i] for i,j in zip(Xpix,Ypix)]

Dist=np.arange(0,np.shape(TransCoord)[0],1)
Temp_trans = [Tz_gr[j, i] for i, j in zip(Xpix, Ypix)]

j=0
for tz in Temp_trans:
    Ts[j]=tz[0]+273.15
    j=j+1
Ttrend=np.array(S_trans)*6e-3
print(Ttrend)
Ts=Ts-Ttrend
Tsnorm=(Ts-min(Ts))/(max(Ts)-min(Ts))
Hnorm=1-(H_trans-min(H_trans))/(max(H_trans)-min(H_trans))



plt.plot(Dist, Tsnorm,'r')
plt.plot(Dist, Hnorm,'b')

plt.xlabel("Distance from Dome C (km)")
plt.ylabel("Brightness temperature (K)")
plt.xlim(0,700)
plt.show()