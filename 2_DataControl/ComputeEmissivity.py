from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']
Tb=Tb[0]

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_MappedonSMOS.nc')
H = np.array(GRISLI.variables['H'])
S = GRISLI.variables['S']
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
TsRACMO=Tz_gr[:,:,0]

TsCroc= netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc')
TsCrocus=TsCroc.variables['TsCrocus'][::-1,:]-273.15

TatDepth=np.zeros(np.shape(TsRACMO))
Tmoy=np.zeros(np.shape(TsRACMO))
Depth=570#570#Ice thickness contributing to the brightness temperature
Layerthick=10 #Layer thickness for interpolation

print("Compute Tmoy")
for i in np.arange(0,np.shape(H)[0]):
    for j in np.arange(0,np.shape(H)[1]):
        if np.shape(np.arange(0,H[i,j],Layerthick))[0]>(Depth/Layerthick):
            #TatDepth[i,j]=np.interp(np.arange(0,H[i,j],Layerthick), (1-Zeta[i,j])*H[i,j], Tz_gr[i,j,:])[Depth//Layerthick]
            Tmoy[i,j]=273.15+np.mean(np.interp(np.arange(0,H[i,j],Layerthick), (1-Zeta[i,j])*H[i,j], Tz_gr[i,j,:])[0:Depth//Layerthick])
        else:
            Tmoy[i, j]=-9999
            #TatDepth[i,j]=-9999

DeltaT=Tmoy-TsRACMO-273.15
TbTs=Tb-TsCrocus-273.15
TmoyTb=Tmoy-Tb

#Calculate dependence of Ts-Tb on Ts-TatDepth
Tb1D=np.reshape(Tb,(np.size(DeltaT),1))
Tmoy1D=np.reshape(Tmoy,(np.size(DeltaT),1))
DeltaT1D=np.reshape(DeltaT,(np.size(DeltaT),1))
TbTs1D=np.reshape(TbTs,(np.size(Tb),1))
TmoyTb1D=np.reshape(TmoyTb,(np.size(Tb),1))

Ts1D=np.reshape(TsCrocus,(np.size(Tb),1))
Tb1D=np.reshape(Tb,(np.size(Tb),1))
Mask1D=np.reshape(Mask,(np.size(Tb),1))

DeltaT1D_=DeltaT1D[DeltaT1D>-50]
TbTs1D_=TbTs1D[DeltaT1D>-50]
TmoyTb1D_=TmoyTb1D[DeltaT1D>-50]
Tmoy1D_=Tmoy1D[DeltaT1D>-50]
TbTs1D__=TbTs1D_[TmoyTb1D_<220]
DeltaT1D__=DeltaT1D_[TmoyTb1D_<220]
TmoyTb1D__=TmoyTb1D_[TmoyTb1D_<220]
Tmoy1D__=Tmoy1D_[TmoyTb1D_<220]

#Least square resolution
# We're solving Ax = B
A = np.column_stack([np.ones(len(DeltaT1D__)), DeltaT1D__])
B = TbTs1D__

# Solve the system of equations.
result,_,_,_= np.linalg.lstsq(A,B, rcond=None)
a, b = result
print(result)

#Correction of the regression model so that consecutive emissivity <= 1
a=0

#2 ways of computing emissivity

#1) With regression model
TbTsTrend1D=a+b*(DeltaT1D)
TbTsModul1D=TbTs1D-TbTsTrend1D
#Reshape data, and correct for the difference between the Ts datasets
TbTsModul=np.reshape(TbTsModul1D,np.shape(TatDepth))+b*(TsCrocus-TsRACMO)
TbTsTrend=np.reshape(TbTsTrend1D,np.shape(TatDepth))-b*(TsCrocus-TsRACMO)
Emissivity=1+TbTsModul/(TbTsTrend/b+TsCrocus)

#2) From the depth where regression coefficient = 1
TbTsModul=np.reshape(TbTsModul1D,np.shape(TatDepth))+b*(TsCrocus-TsRACMO)
TbTsTrend=DeltaT
TbTsModul=TbTs-TbTsTrend-b*(TsCrocus-TsRACMO)
Emissivity=Tb/(Tmoy+b*(TsCrocus-TsRACMO))

Em1D=np.reshape(Emissivity, (np.size(Tb),1))

TbTsModul=TbTsModul*(4-np.array(Mask))/3
Emissivity=Emissivity*(4-np.array(Mask))/3
TbTsTrend=TbTsTrend*(4-np.array(Mask))/3

#Scatter
cmap = plt.cm.get_cmap("coolwarm")
norm=mpl.colors.Normalize(vmin=0.95, vmax=1.05)
plt.scatter(DeltaT1D,TbTs1D, s=5, c=Em1D, cmap=cmap, norm=norm, edgecolors='none')
plt.plot([-10,20],[a+b*(-10),a+20*b],color='r', lw=0.5)
plt.xlabel("Tmoy - Ts (K)")
plt.ylabel("Tb - Ts (K)")
plt.xlim(-10,10)
plt.ylim(-20,20)
plt.colorbar()
plt.grid()
plt.show()

# Geographic plot
fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)
cmap = mpl.cm.spectral
norm = mpl.colors.Normalize(vmin=-15, vmax=15)
myplot=ax[0].pcolormesh(TbTs, cmap=cmap, norm=norm)
ax[0].set_title('Tb-Ts')
ax[1].pcolormesh(TbTsTrend, cmap=cmap, norm=norm)
ax[1].set_title('Regression : Tmoy-Ts')
ax[2].pcolormesh(TbTsModul, cmap=cmap, norm=norm)
ax[2].set_title('Obs-Reg : Tb-Tmoy')
plt.axis('equal')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
cbar = fig.colorbar(myplot, cax=cbar_ax, ticks=np.arange(-15, 16, 5))
cbar.set_label('Kelvin difference (K)', rotation=270)
cbar.ax.set_xticklabels(['-15', '-5', '5', '15'])
plt.show()

fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0.95, vmax=1)
ax.pcolormesh(Emissivity, cmap=cmap, norm=norm)
plt.show()

#Write output NetCDF file
outfile = r'../../SourceData/WorkingFiles/Emissivity.nc'
cols = np.shape(S)[1]
rows = np.shape(S)[0]
dsout = netCDF4.Dataset(outfile, 'w', clobber=True)
Yout = dsout.createDimension('y', rows)
Yout = dsout.createVariable('y', 'f4', ('y',))
Yout.standard_name = 'y'
Yout.units = 'm'
Yout.axis = "Y"
Yout[:] = Y[:,0]+25000
Xout = dsout.createDimension('x', cols)
Xout = dsout.createVariable('x', 'f4', ('x',))
Xout.standard_name = 'x'
Xout.units = 'm'
Xout.axis = "X"
Xout[:] = X[0,:]-25000

dsout.createVariable('Croc-RACMO','float64',('y','x'))
dsout.variables['Croc-RACMO'][:] = np.array(TsCrocus[::-1,:]-TsRACMO[::-1,:])
dsout.createVariable('Ts-Tb','float64',('y','x'))
dsout.variables['Ts-Tb'][:] = np.array(TbTs[::-1,:])
dsout.createVariable('TbTsTrend','float64',('y','x'))
dsout.variables['TbTsTrend'][:] = np.array(TbTsTrend[::-1,:])
dsout.createVariable('TbTsModul','float64',('y','x'))
dsout.variables['TbTsModul'][:] = np.array(TbTsModul[::-1,:])
dsout.createVariable('Emissivity','float64',('y','x'))
dsout.variables['Emissivity'][:] = np.array(Emissivity[::-1,:])
crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
