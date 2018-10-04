#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import netCDF4
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import NC_Resources as ncr
from osgeo import gdal
#np.set_printoptions(threshold=np.nan)

#Control data with a plot - Data is a 2D field
def PlotData(Data, min, max):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    norm = mpl.colors.Normalize(vmin=min, vmax=max)
    cmap = mpl.cm.spectral
    myplot = ax.pcolormesh(Data, cmap=cmap, norm=norm)
    cbar = fig.colorbar(myplot, ticks=np.arange(min, max, (max-min)/5))
    cbar.set_label('Emissivity', rotation=270)
    cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
    plt.show()

#Correct the temperature considering new ice thickness data
def Correct_Tz(H, HBedmap, TsGR, TsCrocus, Tbasal, Tbasalcorr, Tz,GHF, DefHeat, Nx, Ny, Nz):

    Mask = np.zeros(np.shape(H))
    Mask[H >= 10] = 1.0

    # Compute the ratio
    Rapp = np.zeros(np.shape(H))
    Rapp[H > 0] = (Tbasal[H > 0] - TsGR[H > 0]) / (H[H > 0] * GHF[H > 0])

    # Normalize T(z)
    print("Normalization")
    ReducedT = np.zeros(np.shape(Tz))
    for i in np.arange(0, Nx - 1, 1):
        for j in np.arange(0, Ny - 1, 1):
            ReducedT[:, i, j] = (Tz[:, i, j] - Tz[Nz - 1, i, j]) / (Tz[0, i, j] - Tz[Nz - 1, i, j])

    # Show an exemple here:
    [[plt.plot(ReducedT[:, i, j], 1 / (Nz - 1) * np.arange(Nz - 1, -1, -1), linewidth=0.1, c='k') for i in
      np.arange(38, 41, 1)] for j in np.arange(110, 113, 1)]
    plt.show()

    # Critical ice thickness:
    # Considering that Rapp=3e-4 where basal ice is cold
    TypicalRapp = 3e-4
    Hc = np.zeros(np.shape(H))
    Hc[H > 0] = (Tbasal[H > 0] - TsCrocus[H > 0]) / TypicalRapp / (GHF[H > 0] + DefHeat[H > 0])
    Hc[Tbasalcorr == 0] = 0

    # Compute a new basal temperature, depending on the new H (Bedmap2)
    print("New T basal")
    NewTbasal = np.zeros(np.shape(H))
    for i in np.arange(0, Nx - 1, 1):
        for j in np.arange(0, Ny - 1, 1):
            if H[i, j] < Hc[i, j]:
                NewTbasal[i, j] = max(-273.15, Rapp[i, j] * HBedmap[i, j] * GHF[i, j] / TsCrocus[
                    i, j])  # Avoid dummy temperatures
            else:
                NewTbasal[i, j] = Tbasal[i, j]

    # Here assign the average local shape profile to each point
    print("Assign local average shape")
    MeanReducedT = np.zeros(np.shape(Tz))
    for i in np.arange(2, Nx - 3, 1):
        for j in np.arange(2, Ny - 3, 1):
            if np.sum(Mask[i - 1:i + 1, j - 1:j + 1]) != 0:  # Prevent from considering dummy ocean profiles
                MeanReducedT[:, i, j] = np.sum(ReducedT[:, i - 1:i + 1, j - 1:j + 1] * Mask[i - 1:i + 1, j - 1:j + 1],
                                               axis=(1, 2)) / np.sum(Mask[i - 1:i + 1, j - 1:j + 1])

    # Rebuild an updated T(z) from the new local shape
    print("Build the new T profile")
    NewT = np.zeros(np.shape(Tz))
    for i in np.arange(0, Nx-1, 1):
        for j in np.arange(0, Ny-1, 1):
            NewT[:, i, j] = MeanReducedT[:, i, j] * (TsCrocus[i, j] - NewTbasal[i, j]) + NewTbasal[i, j]

    return NewT

#Import GRISLI ice temperature
TGRISLI = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/T3D-AN40C006-k0.nc')
Tz = TGRISLI.variables['T']
TsGR=Tz[0]
Tbasal=Tz[20]

# Import GRISLI velocities and other variables
Model = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/AN40C006_class01_present.nc')
Nx = Model.dimensions['y'].size
Ny = Model.dimensions['x'].size
Nz = Model.dimensions['z'].size
Tbasalcorr = Model.variables['Tb'][0] #Basal T corrected by the melting T
H = Model.variables['H'][0]
GHF=Model.variables['ghf'][0]
Uxbar=Model.variables['Uxbar_maj'][0]
Uybar=Model.variables['Uybar_maj'][0]
Ubar=(Uxbar**2+Uybar**2)**0.5
Taub=1e-3*Model.variables['phid'][0]
DefHeat=Ubar*Taub #Heat caused by ice deformation

nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(Model)

# Import Bedmap2 ice thickness
dataset = gdal.Open('../../SourceData/GRISLI/HBedmap_GRISLI.tif', gdal.GA_ReadOnly)
for x in range(1, dataset.RasterCount + 1):
    band = dataset.GetRasterBand(x)
    HBedmap = band.ReadAsArray()[::-1,:]
HBedmap=HBedmap[::-1,:]

#Import surface temperature Crocus
TSurf= netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/TsCrocus_40km.nc')
TsCrocus=TSurf.variables['TsCrocus'][:,:]

NewT=Correct_Tz(H,HBedmap, TsGR, TsCrocus, Tbasal, Tbasalcorr, Tz,GHF, DefHeat, Nx, Ny, Nz)

#Add normalized depth information
Zeta1=[np.linspace(0,1,21) for _ in range(Ny)]
Zeta2=[np.stack(Zeta1, axis=-1) for _ in range(Nx)]
Zeta=np.stack(Zeta2, axis=-1)

#Create new NetCDF fil for updated temperature data
nc_out = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/Corrected_Tz_GRISLI.nc', 'w', format='NETCDF4')
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