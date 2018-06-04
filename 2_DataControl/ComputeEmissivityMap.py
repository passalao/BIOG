from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

def ComputeEmissivity(Zeta, H, Tz_gr, Tb, Ts):
    Teff=np.zeros(np.shape(H))
    Emissivity=np.zeros(np.shape(H))

    for i in np.arange(0, np.shape(H)[0]):
        for j in np.arange(0, np.shape(H)[1]):
            if Ts[i, j] > -55 and Ts[i, j] < -50 :
                Depth=1000
            if Ts[i, j] > -50 and Ts[i, j] < -45:
                Depth = 600
            if Ts[i, j] > -45 and Ts[i, j] < -40 :
                Depth=375
            if Ts[i, j] > -40 and Ts[i, j] < -35:
                Depth = 325
            if Ts[i, j] > -35:
                Depth = 250
            Teff[i,j]=273.15 + sum(Tz_gr[i, j] * np.exp(-(1 - Zeta[i, j]) * H[i, j] / Depth) * 0.05 * H[i, j] / Depth)/ sum(np.exp(-(1 - Zeta[i, j]) * H[i, j] / Depth) * 0.05 * H[i, j] / Depth)
            Emissivity[i,j]=Tb[i,j]/Teff[i,j]
    return Emissivity

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']
Tb=Tb[0]

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4Ts.nc')
H = np.array(GRISLI.variables['H'])
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

Emissivity=ComputeEmissivity(Zeta, H, Tz_gr, Tb, Ts)

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
