from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import NC_Resources as ncr
import datetime as dt  # Python standard library datetime  module

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_Antarctica_MIR_CDF3TAD_20100112-20180531_52.5deg_epsg6932_nearest.nc')
Tb = Obs.variables['BT_V']
X = Obs.variables['x_spstere']
Y = Obs.variables['y_spstere']
ncols=np.size(X[0,:])
nrows=np.size(Y[:,0])

nc_obsattrs, nc_obsdims, nc_obsvars = ncr.ncdump(Obs)

time=Obs.variables['time'][:]
RefTime=dt.date(2009,1,1)
EndTime=dt.date(2018,5,31)

MonthlyAveragedTbV=np.zeros((12,nrows,ncols))#Dimensions : 12 months, nx, ny, and a field to store the number of useful dates
Lowerbound=min(np.reshape(Tb,(1,np.size(Tb))))[0]
print(Lowerbound)

NbUsefulDates=np.zeros((12,nrows,ncols))

Offset=376
i=Offset
for t in np.arange(376,3437,1):#3437
    Date=RefTime + dt.timedelta(days=i)
    print(Date)
    MonthlyAveragedTbV[Date.month-1,:,:][Tb[i-Offset,:,:]!=Lowerbound]=MonthlyAveragedTbV[Date.month-1,:,:][Tb[i-Offset,:,:]!=Lowerbound]+Tb[i-Offset,:,:][Tb[i-Offset,:,:]!=Lowerbound]
    Mask=Tb[i-Offset,:,:]-Lowerbound
    Mask=np.array([[min(Mask[i,j],1) for j in np.arange(0,ncols,1)] for i in np.arange(0,nrows,1)])
    NbUsefulDates[Date.month-1,:,:]=NbUsefulDates[Date.month-1,:,:]+Mask
    i=i+1

print("Summation finished")

for m in np.arange(0,12,1):
    MonthlyAveragedTbV[m,:,:]=MonthlyAveragedTbV[m,:,:]/NbUsefulDates[m,:,:]
    #MonthlyAveragedTbV[m,:,:][NbUsefulDates[m,:,:]!=0]=MonthlyAveragedTbV[m,:,:][NbUsefulDates[m,:,:]!=0]/NbUsefulDates[m,:,:][NbUsefulDates[m,:,:]!=0]
    #MonthlyAveragedTbV[m,:,:][NbUsefulDates[m,:,:]==0]=9999.0
MonthlyAveragedTbV=np.reshape(MonthlyAveragedTbV, (1,np.size(MonthlyAveragedTbV)))[0]
print(MonthlyAveragedTbV)
MonthlyAveragedTbV[np.isnan(MonthlyAveragedTbV)]=0
MonthlyAveragedTbV=np.reshape(MonthlyAveragedTbV, np.shape(NbUsefulDates))

Sigmas=[[np.std(MonthlyAveragedTbV[:,i,j]) for j in np.arange(0,224,1)] for i in np.arange(0,200,1)]
print(Sigmas)

# Geographic plot
fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0, vmax=1)
cmap = mpl.cm.spectral
myplot=ax.pcolormesh(Sigmas, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0, 1.01, 0.5))
cbar.ax.set_xticklabels(['0', '0.5','1'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.savefig("../../SourceData/SMOS/SMOS_StDev_Monthly.png")
plt.show()

#norm = mpl.colors.Normalize(vmin=200, vmax=0)
#plt.pcolormesh(MonthlyAveragedTbV[1,:,:], cmap=cmap, norm=norm)

#Export NectCDF file
nc_averaged = Dataset('../../SourceData/WorkingFiles/TimeAveragedSMOS.nc', 'w', format='NETCDF4')
nc_averaged.description = "SMOS data averaged over 8 years, month by month"

#Create NetCDF dimensions
nc_averaged.createDimension('m', 12)
nc_averaged.createDimension('x', 200)
nc_averaged.createDimension('y', 224)

nc_averaged.createVariable('TbV', 'float64', ('m','x','y'))
nc_averaged.variables['TbV'][:,:,:] = MonthlyAveragedTbV[:,:,:]