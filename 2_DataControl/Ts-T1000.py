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
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb_Obs = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_MappedonSMOS.nc')
H = np.array(GRISLI.variables['H'])
S = GRISLI.variables['S']
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
TsRACMO=Tz_gr[:,:,0]

# Import temperature data
GRISLIini = netCDF4.Dataset('../../SourceData/WorkingFiles_GRISLIini/GRISLIMappedonSMOS.nc')
Tz_grini = GRISLIini.variables['T']
TsComiso=Tz_grini[:,:,0]

TbTs= netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc')
TsCrocus=TbTs.variables['TsCrocus'][::-1,:]

Tb_Obs=Tb_Obs[0]

GapTsRACMOTb=TsRACMO-Tb_Obs+273.15
GapTsRACMOTb=GapTsRACMOTb*(4-np.array(Mask))/3
GapTsComisoTb=TsComiso-Tb_Obs+273.15
GapTsComisoTb=GapTsComisoTb*(4-np.array(Mask))/3
GapTsCrocusTb=TsCrocus-Tb_Obs
GapTsCrocusTb=GapTsCrocusTb*(4-np.array(Mask))/3

T1000=np.zeros(np.shape(TsCrocus))
T500=np.zeros(np.shape(TsCrocus))

print("Compute T1000")
for i in np.arange(0,np.shape(H)[0]):
    for j in np.arange(0,np.shape(H)[1]):
        if np.shape(np.arange(0,H[i,j],100))[0]>10:
            T1000[i,j]=np.interp(np.arange(0,H[i,j],100), (1-Zeta[i,j])*H[i,j], Tz_gr[i,j,:])[10]
        else:
            T1000[i,j]=-9999
DeltaT1000RACMO=TsRACMO-T1000
DeltaT1000Crocus=TsCrocus-T1000

print("Compute T500")
for i in np.arange(0,np.shape(H)[0]):
    for j in np.arange(0,np.shape(H)[1]):
        if np.shape(np.arange(0,H[i,j],100))[0]>5:
            T500[i,j]=np.interp(np.arange(0,H[i,j],100), (1-Zeta[i,j])*H[i,j], Tz_gr[i,j,:])[5]
        else:
            T500[i,j]=-9999
DeltaT500RACMO = TsRACMO-T500
DeltaT500Crocus=TsCrocus-273.15-T500

#Calculate dependence of Tb on T1000
DeltaT1000_1D=np.reshape(DeltaT1000RACMO,(np.size(DeltaT1000RACMO),1))
TsTb1D=np.reshape(TsCrocus-Tb_Obs,(np.size(Tb_Obs),1))
#TsCrocus1D=np.reshape(TsCrocus,(np.size(Tb_Obs),1))

DeltaT1000_1D_=DeltaT1000_1D[DeltaT1000_1D<1000]
TsTb1D_=TsTb1D[DeltaT1000_1D<1000]
#TsCrocus1D_=TsCrocus1D[T1000_1D>-9999]
print(min(DeltaT1000_1D),max(DeltaT1000_1D))


'''#Least square resolution
# We're solving Ax = B
C = np.column_stack([np.ones(len(DeltaT1000_1D_)), DeltaT1000_1D_])# TsCrocus1D_])
D = TsTb1D_
# Solve the system of equations.
result,_,_,_= np.linalg.lstsq(C, D, rcond=-1)
a, b = result
print(result)
TsTbTrend1D=a+b*(DeltaT1000_1D)
TbModul1D=TsTb1D-TsTbTrend1D
TbModul=np.reshape(TbModul1D,np.shape(T1000))
TbTrend=np.reshape(TsTbTrend1D,np.shape(T1000))'''

# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.magma_r
#norm = mpl.colors.Normalize(vmin=-15, vmax=15)
norm = mpl.colors.Normalize(vmin=-300, vmax=300)
myplot = plt.pcolormesh(DeltaT500RACMO, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0, 10,1))
#cbar = fig.colorbar(myplot, ticks=np.arange(210, 270, 8))
cbar.ax.set_xticklabels(['0', '10'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.show()



#Write output NetCDF file
outfile = r'../../SourceData/WorkingFiles/Deltas_T_et_Tb.nc'
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

dsout.createVariable('TsComiso-TbSMOS','float64',('y','x'))
dsout.variables['TsComiso-TbSMOS'][:] = np.array(GapTsComisoTb[::-1,:])
dsout.createVariable('TsRACMO-TbSMOS','float64',('y','x'))
dsout.variables['TsRACMO-TbSMOS'][:] = np.array(GapTsRACMOTb[::-1,:])
dsout.createVariable('TsCrocus-TbSMOS','float64',('y','x'))
dsout.variables['TsCrocus-TbSMOS'][:] = np.array(GapTsCrocusTb[::-1,:])
dsout.createVariable('TsRACMO-T500','float64',('y','x'))
dsout.variables['TsRACMO-T500'][:] = np.array(DeltaT500RACMO[::-1,:])
dsout.createVariable('TsRACMO-T1000','float64',('y','x'))
dsout.variables['TsRACMO-T1000'][:] = np.array(DeltaT1000RACMO[::-1,:])
dsout.createVariable('TsCrocus-T1000','float64',('y','x'))
dsout.variables['TsCrocus-T1000'][:] = np.array(DeltaT1000Crocus[::-1,:])
dsout.createVariable('TsCrocus-T500','float64',('y','x'))
dsout.variables['TsCrocus-T500'][:] = np.array(DeltaT500Crocus[::-1,:])
dsout.createVariable('TbTrend','float64',('y','x'))
dsout.variables['TbTrend'][:] = np.array(TbTrend[::-1,:])
dsout.createVariable('TbModul','float64',('y','x'))
dsout.variables['TbModul'][:] = np.array(TbModul[::-1,:])
crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
