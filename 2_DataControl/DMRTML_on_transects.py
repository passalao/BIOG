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

#Variables
Transect="DC2McM"#DC2McM, DC2VL or DC2Triangle
TempSource="GRISLI" #Robin or GRISLI

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
Tb=[Tb_Obs[0, j, i] for i,j in zip(Xpix,Ypix)]

# Import temperature data
if TempSource=="Robin":
    FileName="transect_"+Transect+"_Tprof_RobinCrocusGhf60.txt"
    Temp_trans=np.genfromtxt("../../SourceData/Macelloni2016/Macelloni2016_DataTransects/"+FileName, dtype=None, delimiter=" ")
    Temp_trans=Temp_trans-273.15
elif TempSource=="GRISLI":
    GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
    H = GRISLI.variables['H']
    Zeta = GRISLI.variables['Zeta']
    B=GRISLI.variables['B']
    Tz_gr = GRISLI.variables['T']
    Temp_trans=[Tz_gr[j, i] for i,j in zip(Xpix,Ypix)]
    H_trans=[H[j, i] for i,j in zip(Xpix,Ypix)]
    B_trans = [B[j, i] for i, j in zip(Xpix, Ypix)]

# Import Surface temperature: for comparison
TempHelene = netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc')
TsH = TempHelene.variables['TsCrocus']
TsH = TsH[::-1,:]
TstransHelene=[TsH[int(j),int(i)] for i,j in zip(Xpix, Ypix)]

Dist=np.arange(0,np.shape(Temp_trans)[0],1)
Tb_modobs=np.zeros(np.shape(Temp_trans)[0])
Ts=np.zeros(np.shape(Temp_trans)[0])

print('DMRTML computation')
j=0

for tz in Temp_trans:
    if TempSource=="Robin":
        imax=np.where(tz == max(tz))[0][0]
        Thick=imax*50
        tz=tz[0:imax]
    elif TempSource == "GRISLI":
        Thick=H_trans[j]
    Ts[j]=tz[0]+273.15
    Tb_modobs[j] = BIOG.fun.GetTb(tz, Thick, BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams,BIOG.var.Perm, BIOG.var.RTModel,0)
    j=j+1


#plt.plot(Dist, TstransHelene, c='k',lw='1.5')
plt.plot(Dist, B_trans,c='k')
plt.plot(Dist, Ts)
plt.plot(Dist, Tb, c="mediumseagreen", lw='1.5')
plt.plot(Dist, Tb_modobs, c="orangered", lw='1.5') #royalblue for Matzler
plt.grid()
plt.xlabel("Distance from Dome C (km)")
plt.ylabel("Brightness temperature (K)")
plt.xlim(0,700)
plt.ylim(210,235)
plt.show()