from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

def ComputeDepth(Zeta, H, Tz_gr, Tb, TsCorr):
    Depth = np.zeros(np.shape(H))
    e=0.968
    for i in np.arange(0, np.shape(H)[0]):
        for j in np.arange(0, np.shape(H)[1]):
            if Tb[i,j]/e-273.15<min(Tz_gr[i,j])+TsCorr[i,j]:
                #emissivity cannot be 1, the point is not concerned
                Depth[i, j] = -9999
                #TODO; compute something here

            else:
                Found=False
                k=0
                Teff=Tz_gr[i,j,0]
                for t in Tz_gr[i,j,:]:
                    Teff=(k+1)*(Teff+TsCorr[i,j])/(k+2)+(t+TsCorr[i,j])/(k+2)-TsCorr[i,j]
                    #if H[i,j]>2000:
                    #    print(t, Tz_gr[i,j,0],Teff+TsCorr[i,j], Tb[i,j]-273.15, k, Found)
                    #    print(Teff+TsCorr[i,j]- Tb[i,j]+273.15)
                    if Tb[i,j]/e-273.15<Teff+TsCorr[i,j] and Found==False:
                        #print(Teff+TsCorr[i,j]-Tb[i,j]+273.15, k)
                        Found=True
                        Depth[i,j]=(1-Zeta[i,j,k])*H[i,j]
                    k=k+1
    return Depth

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

TsCorr=TsRef-TsRACMO
TatDepth=np.zeros(np.shape(TsRACMO))
Tmoy=np.zeros(np.shape(TsRACMO))
Layerthick=10 #Layer thickness for interpolation

print("Compute Depth")
Depth=ComputeDepth(Zeta, H, Tz_gr, Tb, TsCorr)

fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0, vmax=1000)
cmap = mpl.cm.spectral
myplot=ax.pcolormesh(Depth, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0, 1000, 100))
cbar.set_label('Emissivity', rotation=270)
cbar.ax.set_xticklabels(['0','1000'])
plt.show()

#Write output NetCDF file
outfile = r'../../SourceData/WorkingFiles/PeneDepthE968.nc'
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

dsout.createVariable('PeneDepth','float64',('y','x'))
dsout.variables['PeneDepth'][:] = np.array(Depth[::-1,:])
crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
