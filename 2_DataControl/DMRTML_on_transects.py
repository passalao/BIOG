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
Transect="DC2Triangle"#DC2McM, DC2VL or DC2Triangle
TempSource=["GRISLI", "Robin"] #Robin or GRISLI

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_57.5deg_xy.nc')
Tb_Obs = Obs.variables['BT_V']
X=Obs.variables['x_ease2']
Y=Obs.variables['y_ease2']
Lx=25e3
Ly=25e3

#Import transect coordinates
TransCoord=np.loadtxt("../../SourceData/Transects/"+Transect+".csv", comments="#", delimiter=" ")
#TransCoord=np.loadtxt("../../SourceData/Macelloni2016/Macelloni2016_DataTransects/"+Transect+".csv", comments="#", delimiter=" ")
wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj(init="EPSG:6932")
Xs, Ys = pyproj.transform(wgs84, StereoPol, TransCoord[:,0], TransCoord[:,1])
Xpix=(Xs-X[0][0])//Lx+1
Ypix=(Ys-Y[-1][-1])//Ly+1
Tb=[Tb_Obs[0, j, i] for i,j in zip(Xpix,Ypix)]

# Import temperature data
for t in TempSource:
    if t=="Robin":
        continue
        FileName="transect_"+Transect+"_Tprof_RobinCrocusGhf60.txt"
        Temp_trans=np.genfromtxt("../../SourceData/Macelloni2016/Macelloni2016_DataTransects/"+FileName, dtype=None, delimiter=" ")
        Temp_trans=Temp_trans-273.15
    elif t=="GRISLI":
        GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S004_1_MappedonSMOS.nc')
        H = GRISLI.variables['H']
        S = GRISLI.variables['S']
        Zeta = GRISLI.variables['Zeta']
        B=GRISLI.variables['B']
        Tz_gr = GRISLI.variables['T']
        H_trans=[H[j, i] for i,j in zip(Xpix,Ypix)]
        S_trans=[S[j, i] for i,j in zip(Xpix,Ypix)]
        GradS = np.array([(S_trans[j + 2] - S_trans[j - 2])/ 5e4 for j in np.arange(2, np.size(S_trans)-2, 1)])

        #Temp_trans=[Tz_gr[j, i]-6e-3*((s-2450)) for i,j,s in zip(Xpix,Ypix, S_trans)]
        Temp_trans=[Tz_gr[j, i] for i,j in zip(Xpix,Ypix)]
        #H_trans=[H[j, i] for i,j in zip(Xpix,Ypix)]
        B_trans = [B[j, i] for i, j in zip(Xpix, Ypix)]

    # Import Surface temperature: for comparison
    TempHelene = netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc')
    TsH = TempHelene.variables['TsCrocus']
    TsH = TsH[::-1,:]
    TstransHelene=[TsH[int(j),int(i)] for i,j in zip(Xpix, Ypix)]

    # Import emissivity, computed from regression model
    Emissivity = netCDF4.Dataset('../../SourceData/WorkingFiles/Emissivity.nc')
    E = Emissivity.variables['Emissivity']
    E=E[::-1,:]
    Etrans=[E[int(j),int(i)] for i,j in zip(Xpix, Ypix)]
    Dist=np.arange(0,np.shape(TstransHelene)[0],1)
    Tb_modobs=np.zeros(np.shape(TstransHelene)[0])
    Ts=np.zeros(np.shape(TstransHelene)[0])

    #TODO : corriger Dist en fonction des vraies coordonnées

    print('DMRTML computation: ',t)
    j=0

    for tz in Temp_trans:
        if t=="Robin":
            continue
            imax=np.where(tz == max(tz))[0][0]
            Thick=imax*50
            tz=tz[0:imax]
        elif t == "GRISLI":
            Thick=H_trans[j]
        Ts[j]=tz[0]+273.15
        #tz=tz+TstransHelene[j]-tz[0]-273.15#T(z) correction from bias of Ts RACMO
        Tb_modobs[j] = BIOG.fun.GetTb(tz, Thick, BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams,BIOG.var.Perm, BIOG.var.RTModel,0)
        Tb_modobs[j]=Tb_modobs[j]*Etrans[j]
        j=j+1
    #plt.plot(Dist, TstransHelene, color='dodgerblue', lw='1.5', label="Ts Crocus")
    if t == "Robin":
        continue
        plt.plot(Dist, Tb_modobs, c="blue", lw='1.5', label='Tb Robin') #royalblue for Matzler
    elif t=="GRISLI":
        plt.plot(Dist, Tb_modobs, c="orangered", lw='1.5', label='Tb GRISLI')  # royalblue for Matzler
        plt.plot(Dist, Ts, color="coral", lw='1.5',label='Ts RACMO')
        #plt.plot(Dist[2:-2], GradS, color="coral", lw='1.5',label='GradS')
        #plt.plot(Dist, S_trans, color="k", lw='1.5',label='S')

plt.plot(Dist, Tb, color="mediumseagreen", lw='1.5', label='Tb SMOS')
plt.grid()
plt.legend()
plt.xlabel("Distance from Dome C (km)")
plt.ylabel("Brightness temperature (K)")
plt.xlim(0,700)

#plt.ylim(210,235)
plt.show()