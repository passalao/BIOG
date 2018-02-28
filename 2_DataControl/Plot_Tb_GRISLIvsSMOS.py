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
Model = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb_SMOSGrid.nc')
ny_Mod = Model.dimensions['y'].size
nx_Mod = Model.dimensions['x'].size
Tb_Mod = Model.variables['Tb']

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
Lon = Obs.variables['lon']
Lat = Obs.variables['lat']
Mask = Obs.variables['mask']

GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
Acc = GRISLI.variables['H']

#Select data in the continent
Tb_Mod=Tb_Mod/np.array(Mask)

offset=0
Tb_Obs=np.array(Tb_Obs)
Tb_Mod=np.array(Tb_Mod)
Acc=np.array(Acc)
Error=Tb_Obs[0]-(Tb_Mod+offset)

'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.seismic
norm = mpl.colors.Normalize(vmin=-10, vmax=10)
myplot = plt.pcolormesh(Error, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-15, 15, 1))
cbar.ax.set_xticklabels(['-15', '0', '15'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
#plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()'''

Tb_Obs=np.reshape(Tb_Obs,(201*225,1))
Tb_Mod=np.reshape(Tb_Mod,(201*225,1))
Acc=np.reshape(Acc,(201*225,1))
Mask=np.reshape(Mask,(201*225,1))

Acc4Plot=np.zeros(201*225)
i=0
acclim=0
for a in Acc:
    if a>acclim:
        Acc4Plot[i]=acclim
    else:
        Acc4Plot[i]=a
    i=i+1

#print(np.reshape(np.array(df_clean2["Tb_Obs"]), (1, np.size(df_clean2["Tb_Obs"]))))

# scatterplot
#plt.figure(figsize=(6.5,6.5))
myplot=plt.scatter(Tb_Obs, Tb_Mod+offset, c=Acc, s=0.1)
cbar=plt.colorbar()
cbar.set_label('Ice thickness (m)', rotation=270, labelpad=15)
#cbar.set_label('Geothermal flux (mW/m2)', rotation=270, labelpad=15)
plt.plot([0, 270], [0, 270], color='b')
#plt.autoscale(True)
#plt.text(252,222,"offset = "+str(offset)+" K")
plt.grid()
plt.axis("equal")
plt.autoscale(True)
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb SMOS (K)')
plt.ylabel('Tb GRISLI+SMRT (K)')
#plt.savefig("../../OutputData/img/Tb_SMOSvsMod_DMRTML.png")
plt.show()

'''#Jointplot
#sns.set(style="white")
pdTb_Obs = pd.Series(Tb_Obs.T[0], name="Tb_Obs")
pdTb_Mod = pd.Series(Tb_Mod.T[0], name="Tb_Mod")
pdMask = pd.Series(Mask.T[0], name="Mask")
pdAcc = pd.Series(Acc.T[0], name="Acc")
df=pd.concat((pdTb_Obs, pdTb_Mod, pdAcc, pdMask), axis=1)
df_clean=df[df["Mask"] != 4.0]
df_clean2=df_clean[df_clean["Tb_Obs"] > 100.]

# Show the joint distribution using kernel density estimation
sns.set_style('whitegrid')
g = sns.jointplot(df_clean2["Tb_Obs"], df_clean2["Tb_Mod"], kind="kde", size=7, aspect=1, space=0)#kind="kde", size=7, space=0)
#g = sns.kdeplot(df_clean2["Tb_Obs"], df_clean2["Tb_Mod"], cmap="Reds", shade=True, shade_lowest=False,size=7, ratio=1)#kind="kde", size=7, space=0)
plt.plot([0, 270], [0, 270], color='b')
#g.set(ylim=(200,270))
#g.set(xlim=(200,270))
#plt.autoscale(True)
#plt.axis("equal")
#plt.xlim(200, 270)
#plt.ylim(200, 270)
#plt.grid()
plt.show()'''
plt.close()


'''
#Once used to reproject data
Lon = np.reshape(Lon, (nx_Obs * ny_Obs, 1))
Lat = np.reshape(Lat, (nx_Obs * ny_Obs, 1))

# GRISLI coordinates
xG = np.linspace(-2.805e6, 3e6, 387)
yG = np.linspace(2.805e6, -2.805e6, 374)

# SMOS coordinates
import mpl_toolkits.basemap.pyproj as pyproj

wgs84 = pyproj.Proj("+init=EPSG:4326")
StereoPol = pyproj.Proj("+init=EPSG:3031")  EPSG:6932 EASE2 grid#
XX, YY = pyproj.transform(wgs84, StereoPol, Lon, Lat)
X = np.reshape(XX, (nx_Obs, ny_Obs))
Y = np.reshape(YY, (nx_Obs, ny_Obs))
Tb_Obs=Tb_Obs[0,:,:]
f = si.interp2d(X[0, :], Y[:, 0], Tb_Obs, kind='cubic')
Tb_SMOS = f(xG, yG)  # Interpolate on the new grid'''