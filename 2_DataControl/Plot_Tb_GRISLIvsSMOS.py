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

# Import Tb data computed from the model
#Model = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb.nc')

Model1 = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb_SMOSGrid_'+BIOG.var.RTModel+'_'+BIOG.var.Perm+'_Rexp.nc')
ny_Mod1 = Model1.dimensions['y'].size
nx_Mod1 = Model1.dimensions['x'].size
Tb_Mod1 = Model1.variables['Tb']

#Model2 = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI_Tb_SMOSGrid_'+'DMRT-ML'+'_'+'Matzler'+'.nc')
#Tb_Mod2 = Model2.variables['Tb']

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']

GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
Acc = GRISLI.variables['Acc']
T = GRISLI.variables['T']
H = GRISLI.variables['H']
Ts=T[:,:,0]
Tground=T[:,:,20]
rho=917
g=9.81
PMP=273.16-0.074*rho*g*(H[:,:]-30)/1e6-0.024
GapToPMP=Tground-PMP+273.15

offset=0.
Tb_Obs=np.array(Tb_Obs)
Tb_Mod1=np.array(Tb_Mod1)
#Tb_Mod2=np.array(Tb_Mod2)
Acc=np.array(Acc)
Ts=np.array(Ts)
Tb_Mod1=Tb_Mod1*(4-np.array(Mask))/3
Tb_Obs=Tb_Obs*(4-np.array(Mask))/3
GapToPMP=GapToPMP*(4-np.array(Mask))/3
Error=Tb_Mod1-Tb_Obs[0]
Error=Error*(4-np.array(Mask))/3
GapTsTb=Ts-Tb_Obs[0]+273.15
GapTsTb=GapTsTb*(4-np.array(Mask))/3

#Tb_Mod2=Tb_Mod2/np.array(Mask)


Tb_Obs=np.reshape(Tb_Obs,(201*225,1))
Tb_Mod1=np.reshape(Tb_Mod1,(201*225,1))
#Tb_Mod2=np.reshape(Tb_Mod2,(201*225,1))
Acc=np.reshape(Acc,(201*225,1))
Ts=np.reshape(Ts,(201*225,1))
Mask=np.reshape(Mask,(201*225,1))

'''Acc4Plot=np.zeros(201*225)
i=0
acclim=0.4
for a in Acc:
    if a>acclim:
        Acc4Plot[i]=acclim
    else:
        Acc4Plot[i]=a
    i=i+1'''

'''outfile = r'../../SourceData/WorkingFiles/SMOSvsModel.nc'
cols = len(X[0,:])
rows = len(Y[:,0])

dsout = netCDF4.Dataset(outfile, 'w', clobber=True)

Yout = dsout.createDimension('y', rows)
Yout = dsout.createVariable('y', 'f4', ('y',))
Yout.standard_name = 'y'
Yout.units = 'm'
Yout.axis = "Y"
Yout[:] = Y[:,0]

Xout = dsout.createDimension('x', cols)
Xout = dsout.createVariable('x', 'f4', ('x',))
Xout.standard_name = 'x'
Xout.units = 'm'
Xout.axis = "X"
Xout[:] = X[0,:]

dsout.createVariable('TbSMOS','float64',('y','x'))
dsout.variables['TbSMOS'][:] = np.array(Tb_Obs[::-1,::-1])
dsout.createVariable('TsGRISLI','float64',('y','x'))
dsout.variables['TsGRISLI'][:] = np.array(Ts[::-1,:])
dsout.createVariable('DeltaToPMP','float64',('y','x'))
dsout.variables['DeltaToPMP'][:] = np.array(GapToPMP[::-1,:])
dsout.createVariable('TsGRISLI-TbSMOS','float64',('y','x'))
dsout.variables['TsGRISLI-TbSMOS'][:] = np.array(GapTsTb[::-1,:])
dsout.createVariable('TbModel','float64',('y','x'))
dsout.variables['TbModel'][:] = np.array(Tb_Mod1[::-1,:])
dsout.createVariable('TbGRISLI-TbSMOS','float64',('y','x'))
dsout.variables['TbGRISLI-TbSMOS'][:] = np.array(Error[::-1,:])

crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
'''

# scatterplot
#myplot=plt.scatter(Tb_Mod1, Tb_Mod2, c="Red", s=5, label="Mätzler")
myplot=plt.scatter(Tb_Obs, Tb_Mod1, c="Red", s=5, label="Mätzler")
#myplot=plt.scatter(Tb_Obs, Tb_Mod2+offset, c="Darkgreen", s=5, label="Tiuri")
#cbar=plt.colorbar()
#cbar.set_label('Geothermal flux (mW/m2)', rotation=270, labelpad=15)
#cbar.set_label('Accumulation (m/a)', rotation=270, labelpad=15)
plt.plot([0, 270], [0, 270], color='b')
plt.grid()
plt.axis("equal")
plt.autoscale(True)
plt.legend()
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb SMOS (K)')
plt.ylabel('Tb GRISLI+'+BIOG.var.RTModel+'(K)')
#plt.savefig("../../OutputData/img/Tb_SMOSvsMod_"+BIOG.var.RTModel+".png")
plt.show()

'''#plt.figure(figsize=(6.5,6.5))
plt.scatter(Ts+273.15, Tb_Mod1, c="Red", s=10)
plt.scatter(Ts+273.15, Tb_Mod2, c="Green", s=10)
plt.plot([0, 270], [0, 270], color='b')
plt.grid()
plt.axis("equal")
plt.autoscale(True)
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb Mätzler (K)')
plt.ylabel('Tb Tiuri (K)')
#plt.savefig("../../OutputData/img/Tb_SMOSvsMod_DMRTML.png")
plt.show()'''

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

'''# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.magma_r
norm = mpl.colors.Normalize(vmin=0, vmax=20)
myplot = plt.pcolormesh(Error, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0, 21, 5))
cbar.ax.set_xticklabels(['-15', '0', '15'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
#plt.savefig("../../OutputData/img/Error_SMOS-Mod_DMRTML.png")
plt.show()'''

