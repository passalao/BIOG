#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
from netCDF4 import Dataset
import netCDF4, BIOG
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
import scipy.interpolate as si
import NC_Resources as ncr
from osgeo import gdal
np.set_printoptions(threshold=np.nan)

def PlotData(Data, min, max):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    norm = mpl.colors.Normalize(vmin=min, vmax=max)
    cmap = mpl.cm.spectral
    myplot = ax.pcolormesh(Data, cmap=cmap, norm=norm)
    cbar = fig.colorbar(myplot, ticks=np.arange(min, max, (max-min)/5))
    cbar.set_label('Emissivity', rotation=270)
    cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
    plt.savefig("../../OutputData/img/InvertingEmissDepth/Emissivity_DescentGrad_nDim.png")
    plt.show()

# Import GRISLI data
Model = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/AN40C006_class01_present.nc')
ny_Mod = Model.dimensions['y'].size
nx_Mod = Model.dimensions['x'].size
nz_Mod = Model.dimensions['z'].size
nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(Model)
Tbasalcorr = Model.variables['Tb'][0]
H = Model.variables['H'][0]
GHF=Model.variables['ghf'][0]
Uxbar=Model.variables['Uxbar_maj'][0]
Uybar=Model.variables['Uybar_maj'][0]
Ubar=Uxbar**2+Uybar**2
Taub=1e-3*Model.variables['phid'][0]
DefHeat=Ubar*Taub

# Import Bedmap2 ice thickness
dataset = gdal.Open('../../SourceData/GRISLI/HBedmap_GRISLI.tif', gdal.GA_ReadOnly)
for x in range(1, dataset.RasterCount + 1):
    band = dataset.GetRasterBand(x)
    HBedmap = band.ReadAsArray()

#Import GRISLI ice temperature
TGRISLI = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/T3D-AN40C006-k0.nc')
Tz = TGRISLI.variables['T']
TsGR=Tz[0]
Tbasal=Tz[20]
nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(Model)

#Import surface temperature Crocus
TSurf= netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/TsCrocus_40km.nc')
TsCrocus=TSurf.variables['TsCrocus'][:,:]

Rapp=np.zeros(np.shape(H))
Rapp[H>0]=(Tbasal[H>0]-TsGR[H>0])/(H[H>0]*GHF[H>0])

#Calul de la forme de T(z) normalisée
print("Normalization")
ReducedT=np.zeros(np.shape(Tz))
for i in np.arange(0,140,1):
    for j in np.arange(0, 140, 1):
        ReducedT[:,i,j]=(Tz[:,i,j]-Tz[20,i,j])/(Tz[0,i,j]-Tz[20,i,j])

[[plt.plot(ReducedT[:, i, j], 0.05*np.arange(20,-1,-1), linewidth=0.1, c='k') for i in np.arange(38,41,1)] for j in np.arange(110,113,1)]
plt.show()

#Critical ice thickness:
#Considering that Rapp=3e-4 where basal ice is cold
Hc=np.zeros(np.shape(H))
Hc[H>0]=(Tbasal[H>0]-TsCrocus[H>0])/3e-4/(GHF[H>0]+DefHeat[H>0])
Hc[Tbasalcorr==0]=0

#Ici calculer un nouveau Tbasal, en fonction du H de bedmap
print("New T basal")
NewTbasal=np.zeros(np.shape(H))
for i in np.arange(0,140,1):
    for j in np.arange(0, 140, 1):
        if H[i,j]<Hc[i,j]:
            NewTbasal[i,j]=max(-273.15,Rapp[i,j]*HBedmap[i,j]*GHF[i,j]/TsCrocus[i,j])#Avoid dummy temperatures
        else:
            NewTbasal[i,j]=Tbasal[i,j]

#Ici autour de chaque point assigner la forme moyenne des points autour
print("Assign local average shape")
MeanReducedT=np.zeros(np.shape(Tz))
for i in np.arange(2,138,1):
    for j in np.arange(2, 138, 1):
        MeanReducedT[:,i,j]=np.mean(ReducedT[:,i-1:i+1,j-1:j+1], axis=(1,2))

#Reconstruire un T(z) à partir de la forme
print("Build the new T profile")
NewT=np.zeros(np.shape(Tz))
Zeta=np.zeros(np.shape(Tz))
for i in np.arange(0,140,1):
    for j in np.arange(0, 140, 1):
        NewT[:,i,j]=MeanReducedT[:,i,j]*(TsCrocus[i,j]-NewTbasal[i,j])+NewTbasal[i,j]
        Zeta[:,i,j]=0.05*np.arange(0,21,1)

#Create new NetCDF fil for interpolated data
nc_out = Dataset('../../SourceData/GRISLI/Avec_FoxMaule/Corrected_Tz_GRISLI.nc', 'w', format='NETCDF4')
nc_out.description = "GRISLI data mapped on SMOS grid corrected for the ice thickness"

#Create NetCDF dimensions
for dim in nc_moddims:
    dim_name=dim
    if dim=="rows":
        dim_name="x"
    if dim=="cols":
        dim_name="y"
    nc_out.createDimension(dim_name, Model.dimensions[dim].size)

nc_out.createVariable('T', 'float64', ('z','y','x'))
nc_out.variables['T'][:] = NewT
nc_out.createVariable('H', 'float64', ('y','x'))
nc_out.variables['H'][:] = HBedmap [::-1,:]
nc_out.createVariable('Zeta', 'float64', ('z','y','x'))
nc_out.variables['Zeta'][:] = 1-Zeta [:, ::-1,:]
nc_out.close()