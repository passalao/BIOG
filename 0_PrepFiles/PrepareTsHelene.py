#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import NC_Resources as ncr
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from osgeo import gdal

dataset = gdal.Open('../../SourceData/Crocus/TS_CROCUS_Helene.img', gdal.GA_ReadOnly)

print("Driver: {}/{}".format(dataset.GetDriver().ShortName,
                             dataset.GetDriver().LongName))
print("Size is {} x {} x {}".format(dataset.RasterXSize,
                                    dataset.RasterYSize,
                                    dataset.RasterCount))
print("Projection is {}".format(dataset.GetProjection()))

geotransform = dataset.GetGeoTransform()
if geotransform:
    print("Origin = ({}, {})".format(geotransform[0], geotransform[3]))
    print("Pixel Size = ({}, {})".format(geotransform[1], geotransform[5]))

band = dataset.GetRasterBand(1)
print("Band Type={}".format(gdal.GetDataTypeName(band.DataType)))

min = band.GetMinimum()
max = band.GetMaximum()
if not min or not max:
    (min, max) = band.ComputeRasterMinMax(True)
print("Min={:.3f}, Max={:.3f}".format(min, max))

if band.GetOverviewCount() > 0:
    print("Band has {} overviews".format(band.GetOverviewCount()))

if band.GetRasterColorTable():
    print("Band has a color table with {} entries".format(band.GetRasterColorTable().GetCount()))

scanline = band.ReadRaster(xoff=0, yoff=0,
                           xsize=band.XSize, ysize=band.YSize,
                           buf_xsize=band.XSize, buf_ysize=band.YSize,
                           buf_type=gdal.GDT_Float32)
import struct
tuple_of_floats = struct.unpack('f' * band.XSize*band.YSize, scanline)
Ts=np.reshape(np.array(tuple_of_floats),(band.YSize,band.XSize))


'''print("Check origin")
for j in np.arange(0,band.YSize):
    for i in np.arange(0,band.XSize):
        Xp=geotransform[0]+j*geotransform[1]
        Yp=geotransform[3]+i*geotransform[5]
        if Xp==0.0 and Yp==0.0:
            Ts[i,j]=270.
            print("South pole :", i, j)'''

Ts=Ts[::-1,:]
Ts[Ts>270]=-999999.0
Ts[Ts==-999999.0]=0

#Correction for the number of pixels
Corr_NbPixels=np.array(-999999.0*np.zeros((3,band.XSize)))
Ts=np.concatenate((Corr_NbPixels, Ts),axis=0)
Ts=Ts[:-3,6:]

#Import GRISLI data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S004_1_MappedonSMOS.nc')#GRISLIMappedonSMOS.nc')
Acc = GRISLI.variables['Acc']
T = GRISLI.variables['T']
H = GRISLI.variables['H']
UXBAR = GRISLI.variables['UXBAR']
UYBAR = GRISLI.variables['UYBAR']
Ts_gr=T[:,:,0]+273.15
Us=(np.array(UXBAR)**2+np.array(UYBAR)**2)**0.5*4/3
# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']
Tb_Obs=np.array(Tb_Obs)
Tb_Obs=Tb_Obs*(4-np.array(Mask))/3
Tb_Obs=Tb_Obs[0]
GapTsTb=Ts-Tb_Obs
GapTsTb=GapTsTb*(4-np.array(Mask))/3

Ts_gr=Ts_gr*(4-np.array(Mask))/3
Us=Us*(3.9999-np.array(Mask))/3

'''i=79#80#143#86#65
j=93#94#147#165#167
print(Ts[i,j], Tb_Obs[0,i,j])'''
'''
# scatterplot
myplot=plt.scatter(Tb_Obs, Ts, c="Red", s=0.01)
plt.plot([0, 270], [0, 270], color='b')
plt.grid()
plt.axis("equal")
#plt.autoscale(True)
plt.legend()
plt.xlim(200, 270)
plt.ylim(200, 270)
plt.xlabel('Tb SMOS (K)')#Ts Comiso (K)')
plt.ylabel('Ts Crocus (K)')
#plt.savefig("../../OutputData/img/Tb_SMOSvsMod_"+BIOG.var.RTModel+".png")
#plt.show()

# Geographic plot
fig, ax = plt.subplots()
cmap = mpl.cm.magma_r
#norm = mpl.colors.Normalize(vmin=-15, vmax=15)
norm = mpl.colors.Normalize(vmin=210, vmax=270)
#myplot = plt.pcolormesh(GapTsTb, cmap=cmap, norm=norm)
myplot = plt.pcolormesh(Ts, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-15, 15.1,5))
#cbar = fig.colorbar(myplot, ticks=np.arange(210, 270, 8))
cbar.ax.set_xticklabels(['210', '270'])  # vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
#plt.savefig("../../OutputData/img/Error_SMOS-Mod_DMRTML.png")
'''

#Write output NetCDF file
outfile = r'../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc'
cols = len(X[0,:])
rows = len(Y[:,0])
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

dsout.createVariable('TbSMOS','float64',('y','x'))
dsout.variables['TbSMOS'][:] = np.array(Tb_Obs[::-1,:])
dsout.createVariable('TsCrocus','float64',('y','x'))
dsout.variables['TsCrocus'][:] = np.array(Ts[::-1,:])
dsout.createVariable('TsCrocus-TbSMOS','float64',('y','x'))
dsout.variables['TsCrocus-TbSMOS'][:] = np.array(GapTsTb[::-1,:])
dsout.createVariable('TsCrocus-TsComiso','float64',('y','x'))
dsout.variables['TsCrocus-TsComiso'][:] = np.array(Ts[::-1,:]-Ts_gr[::-1,:])
dsout.createVariable('Us','float64',('y','x'))
dsout.variables['Us'][:] = np.array(Us[::-1,:])
crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
