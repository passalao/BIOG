from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

def ComputeTmoy(Zeta, H, Tz_gr, Layerthick, TsRef):
    Tmoy = np.zeros(np.shape(H))
    for i in np.arange(0, np.shape(H)[0]):
        for j in np.arange(0, np.shape(H)[1]):
            #Determine Ts:
            if TsRef[i,j]<-52.5:
                Depth=1350
                #Depth=1600
            if TsRef[i,j]<-50 and TsRef[i,j]>-52.5:
                Depth = 1250
                #Depth = 1250
            if TsRef[i,j]<-47.5 and TsRef[i,j]>-50:
                Depth = 730
                #Depth = 740
            if TsRef[i,j]<-45 and TsRef[i,j]>-47.5:
                Depth = 680
                #Depth = 710
            if TsRef[i,j]<-40 and TsRef[i,j]>-45:
                Depth = 480
                #Depth = 470
            if TsRef[i,j]<-35 and TsRef[i,j]>-40:
                Depth = 160
                #Depth = 580
            if TsRef[i,j]>-35:
                Depth = 40
                #Depth = 270

            if np.shape(np.arange(0, H[i, j], Layerthick))[0] > (Depth / Layerthick):
                Tmoy[i, j] = 273.15 + np.mean(
                    np.interp(np.arange(0, H[i, j], Layerthick), (1 - Zeta[i, j]) * H[i, j], Tz_gr[i, j, :])[
                    0:Depth // Layerthick])
            else:
                Tmoy[i, j] = -9999
    return Tmoy

TsData="Crocus"#"MAR" or "Crocus"

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
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
TsRACMO=Tz_gr[:,:,0]

if TsData=="MAR":
    TsR = netCDF4.Dataset('../../SourceData/WorkingFiles/TsMAR.nc')
    TsRef = TsR.variables['TsMAR'][::-1, :]
elif TsData=="Crocus":
    TsR= netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc')
    TsRef=TsR.variables['TsCrocus'][::-1,:]-273.15

TatDepth=np.zeros(np.shape(TsRACMO))
Tmoy=np.zeros(np.shape(TsRACMO))
Layerthick=10 #Layer thickness for interpolation

print("Compute Tmoy")
Tmoy=ComputeTmoy(Zeta, H, Tz_gr, Layerthick, TsRef)
print(Tmoy)
Emissivity=Tb/(Tmoy+TsRef-TsRACMO)

#Emissivity=Emissivity*(4-np.array(Mask))/3

fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0.95, vmax=1)
cmap = mpl.cm.spectral
myplot=ax.pcolormesh(Emissivity, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0.95, 1.01, 0.01))
cbar.set_label('Emissivity', rotation=270)
cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
plt.show()

#Write output NetCDF file
outfile = r'../../SourceData/WorkingFiles/Emissivity.nc'
cols = np.shape(H)[1]
rows = np.shape(H)[0]
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

dsout.createVariable('Emissivity','float64',('y','x'))
dsout.variables['Emissivity'][:] = np.array(Emissivity[::-1,:])
crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
