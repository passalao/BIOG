from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
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

time=Obs.variables['time'][:]
RefTime=dt.date(2009,1,1)
StartTime=dt.date(2010, 1, 12)
EndTime=dt.date(2018,5,31)

MonthlyAveragedTbV=np.zeros((12,nrows,ncols))#Dimensions : 12 months, nx, ny, and a field to store the number of useful dates
Lowerbound=min(np.reshape(Tb,(1,np.size(Tb))))[0]
NbUsefulDates=np.zeros((12,nrows,ncols))

Offset=376
i=Offset
for t in np.arange(376,400,1):#time:
    Date=RefTime + dt.timedelta(days=i)
    print(Date)

    MonthlyAveragedTbV[Date.month-1,:,:]=MonthlyAveragedTbV[Date.month-1,:,:]+Tb[i-Offset,:,:]
    Mask=Tb[i-Offset,:,:]-Lowerbound
    Mask=[[min(Mask[i,j],1) for j in np.arange(0,ncols,1)] for i in np.arange(0,nrows,1)]
    NbUsefulDates[Date.month-1,:,:]=NbUsefulDates[Date.month-1,:,:]+Mask
    i=i+1

print("Summation finished")

for m in np.arange(0,12,1):
    print(m)
    MonthlyAveragedTbV[m,:,:]=MonthlyAveragedTbV[m,:,:]/NbUsefulDates[m,:,:]

#Export NectCDF file
