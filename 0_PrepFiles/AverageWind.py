from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4
#import matplotlib.pyplot as plt
#import matplotlib as mpl
import numpy as np
#import NC_Resources as ncr
import datetime as dt  # Python standard library datetime  module

# Import SMOS data (for the coordinate, missing in wind source file)
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
NewX=np.arange(X[0,0]-3*25000, X[0,-1]+4*25000, 25000)

# Import Wind data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/Vent/ERAIN-1980-2016-Wind.nc')
Wind = Obs.variables['Wind']

ncols=np.size(Wind[0,0,:])
nrows=np.size(Wind[0,:,0])

 #nc_obsattrs, nc_obsdims, nc_obsvars = ncr.ncdump(Obs)

time=Obs.variables['time'][:]
AveragedWind=np.zeros((nrows,ncols))

'''for i in np.arange(0,nrows,1):
    for j in np.arange(0, ncols, 1):
        print(i,' ', j)
        AveragedWind[i,j]=max(Wind[0:365,i,j])'''
#AveragedWind=np.array([[np.mean(Wind[0:365,i,j]) for j in np.arange(0, ncols,1)] for i in np.arange(0,nrows,1)])
AveragedWind=np.array([[max(Wind[0:365,i,j]) for j in np.arange(0, ncols,1)] for i in np.arange(0,nrows,1)])

'''
fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0, vmax=10)
cmap = mpl.cm.spectral
myplot = ax.pcolormesh(AveragedWind, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0.90, 1.01, 0.02))
cbar.set_label('Emissivity', rotation=270)
cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
plt.show()
'''

#Create NetCDF file
cols = np.shape(AveragedWind)[1]
rows = np.shape(AveragedWind)[0]

outfile = r'../../SourceData/Vent/MaxdWind.nc'
nc_new = netCDF4.Dataset(outfile, 'w', clobber=True)

Yout = nc_new.createDimension('y', rows)
Yout = nc_new.createVariable('y', 'f4', ('y',))
Yout.standard_name = 'y'
Yout.units = 'm'
Yout.axis = "Y"
Yout[:] = Y[:,0]+25000
Xout = nc_new.createDimension('x', cols)
Xout = nc_new.createVariable('x', 'f4', ('x',))
Xout.standard_name = 'x'
Xout.units = 'm'
Xout.axis = "X"
Xout[:] = NewX[:]-25000

nc_new.createVariable("Wind", 'float64', ('y','x'))
nc_new.variables["Wind"][:] = AveragedWind[:,:]
crs = nc_new.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
nc_new.close()