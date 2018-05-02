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

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
Tb_Obs = Obs.variables['BT_V']
Lat=Obs.variables['lat']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Lx=25e3
Ly=25e3
Tb_Obs=Tb_Obs[0]

'''# Import Crocus data
ObsCrocus = netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc')
Ts = ObsCrocus.variables['TsCrocus'][::-1,:]
Ts=Ts-273.15
Ts[Ts==-273.15]=0.0'''

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_MappedonSMOS.nc')
H = np.array(GRISLI.variables['H'])
S = GRISLI.variables['S']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

#Calculate surface gradients
GradS=np.array(np.gradient(S))
Slope=(GradS[0,:,:]**2+GradS[1,:,:]**2)**0.5/5e4
LapSx=np.array(np.gradient(GradS[0,:,:]))
LapSy=np.array(np.gradient(GradS[1,:,:]))
Curv=(LapSx[0,:,:]**2+LapSy[1,:,:]**2)**0.5/5e4
Curv=(LapSx[0,:,:]+LapSy[1,:,:])/2/5e4

#Calculate dependence of Ts on S and Lat
S1D=np.reshape(S,(np.size(S),1))
H1D=np.reshape(H,(np.size(H),1))
Lat1D=np.reshape(Lat,(np.size(Lat),1))
Ts1D=np.reshape(Ts,(np.size(Ts),1))

#Least square resolution
# We're solving Ax = B
A = np.column_stack([np.ones(len(S1D)), S1D, Lat1D])
B = Ts1D
# Solve the system of equations.
result, _, _,_= np.linalg.lstsq(A, B, rcond=None)
a, b, c = result
TsTrend1D=a+b*S1D+c*Lat1D
TsModul1D=Ts1D-TsTrend1D
TsModul=np.reshape(TsModul1D,np.shape(S))
TsTrend=np.reshape(TsTrend1D,np.shape(S))

j=100
i=130
'''
plt.plot(np.arange(1,224,1),1e3*GradS[j,:], label='Grad S')
#plt.plot(np.arange(0,225,1),1e-2*H[j,:], label='H')
plt.plot(np.arange(0,225,1),1e-2*S[j,:], label='S')
plt.plot(np.arange(0,225,1),-Ts[j,:], label='Ts')
plt.plot(np.arange(0,225,1),-TsTrend[j,:], label='Ts trend')
plt.plot(np.arange(0,225,1),TsModul[j,:], label='Ts modul')'''

plt.plot(np.arange(0,201,1),1e4*Curv[:,i], label='Grad S')
#plt.plot(np.arange(0,201,1),1e-2*H[:,i], label='H')
plt.plot(np.arange(0,201,1),1e-2*S[:,i], label='S')
plt.plot(np.arange(0,201,1),-Ts[:,i], label='Ts')
plt.plot(np.arange(0,201,1),-TsTrend[:,i], label='Ts trend')
plt.plot(np.arange(0,201,1),TsModul[:,i], label='Ts modul')
plt.legend()
#plt.show()

#plt.scatter(TsModul[np.array(S)>3000], np.array(H)[np.array(S)>3000], s=0.1)
plt.scatter(TsModul[150:200,70:150], GradS[0,149:199,69:149], s=0.1)
plt.show()

# Geographic plot
fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,figsize=(6.75,6.03))
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=210, vmax=270)
ax[0,0].pcolormesh(Ts+273.15, cmap=cmap, norm=norm)
norm = mpl.colors.Normalize(vmin=210, vmax=270)

myplot2=ax[0,1].pcolormesh(TsTrend+273.15, cmap=cmap, norm=norm)
cmap = mpl.cm.gist_earth_r
norm = mpl.colors.Normalize(vmin=-15, vmax=15)
myplot4=ax[1,0].pcolormesh(TsModul, cmap=cmap, norm=norm)

ax[1,0].set_xlim([0, 225])
ax[1,0].set_ylim([0, 201])
cmap = mpl.cm.gist_earth_r
#norm = mpl.colors.Normalize(vmin=1000, vmax=4000)
norm = mpl.colors.Normalize(vmin=0, vmax=2e-3)
myplot4=ax[1,1].pcolormesh(Slope, cmap=cmap, norm=norm)
#cbar = fig.colorbar(myplot, ticks=np.arange(0, 4000, 500))
#cbar.ax.set_xticklabels(['0', '2000', '4000'])
#ax[1,1].pcolormesh(H, cmap=cmap, norm=norm)
ax[1,1].set_xlim([0, 225])
ax[1,1].set_ylim([0, 201])
# vertically oriented colorbar
#plt.axis('equal')
#plt.autoscale(True)
plt.show()

#Write output NetCDF file
outfile = r'../../SourceData/WorkingFiles/TsModulations_and_H.nc'
cols = np.shape(S)[1]
rows = np.shape(S)[0]
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

dsout.createVariable('Ts','float64',('y','x'))
dsout.variables['Ts'][:] = np.array(Ts[::-1,:])
dsout.createVariable('TsTrend','float64',('y','x'))
dsout.variables['TsTrend'][:] = np.array(TsTrend[::-1,:])
dsout.createVariable('TsModul','float64',('y','x'))
dsout.variables['TsModul'][:] = np.array(TsModul[::-1,:])
dsout.createVariable('Tb','float64',('y','x'))
dsout.createVariable('H','float64',('y','x'))
dsout.variables['H'][:] = np.array(H[::-1,:])
dsout.createVariable('Slope','float64',('y','x'))
dsout.variables['Slope'][:] = np.array(Slope[::-1,:]*1e3)
dsout.createVariable('Curvature','float64',('y','x'))
dsout.variables['Curvature'][:] = np.array(Curv[::-1,:]*1e4)
crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
