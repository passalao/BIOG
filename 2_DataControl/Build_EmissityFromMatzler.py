from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import NC_Resources as ncr
import scipy.optimize as opt
import scipy.interpolate as interp
import scipy.integrate as sciint
np.set_printoptions(threshold=np.nan)

###################################################################################################
#Functions
###################################################################################################
#Computes the effective temperature
def ComputeTeff():
    Freq=1.413e9
    l=299792458/Freq

    Teff=np.zeros(np.shape(H))
    for i in np.arange(0,np.shape(H)[0], 1):
        for j in np.arange(0,np.shape(H)[1], 1):
            if H[i,j]>1.0:
                Perm = Permittivity(Tz_gr[i, j])
                #print(Perm)
                Lz = l / (4 * math.pi * Perm[1, :]) * (Perm[0,:]) ** 0.5
                z = (1 - Zeta[i, j]) * H[i, j]
                dz=0.05*H[i,j]
                #zeta=np.linspace(0,1,21)
                '''IntExp=np.zeros(np.size(Lz))
                for k in np.arange(0,np.size(Lz)-1,1):
                    IntExp[k]=sum(((1/Lz[0:-2]+1/Lz[1:-1])/2*dz)[0:k])'''

                #Computation of effective temperature. It seems OK, a bit low maybe
                #Teff[i,j]=sum((Tz_gr[i,j][0:-2]/Lz[0:-2]*np.exp(-IntExp[0:-2])+(Tz_gr[i,j][1:-1]/Lz[1:-1]*np.exp(-IntExp[1:-1])))/2*dz)

                #Gauss-Legendre integration (default interval is [-1,1])
                deg = 10
                x, w = np.polynomial.legendre.leggauss(deg)
                a=0
                b=1
                t = 0.5 * (x + 1) * (b - a) + a
                Teff[i, j] = sum(w * Integrand(t,Tz_gr[i, j], Lz, H[i, j])) * 0.5 * (b - a)
                print(Teff[i,j])
    return Teff

def Integrand2(zeta2, Lz2,H):
    Linterp2 = interp.interp1d(np.linspace(0, 1, np.size(Lz2)), Lz2)
    L2 = Linterp2(zeta2)
    return H/L2

def Integrand(zeta,Tz,Lz,H):
    Tinterp=interp.interp1d(np.linspace(0,1,np.size(Tz)), Tz)
    Linterp=interp.interp1d(np.linspace(0,1,np.size(Lz)), Lz)
    T=Tinterp(zeta)
    L=Linterp(zeta)

    # Gauss-Legendre integration (default interval is [-1,1])
    deg = 10
    x1, w1 = np.polynomial.legendre.leggauss(deg)
    a1 = 0
    b1 = zeta
    t1 = 0.5 * (x1 + 1) * (b1 - a1) + a1
    IntExp = sum(w1 * Integrand2(t1,L,H) * 0.5 * (b1 - a1))
    #return T *H* np.exp(-zeta*H / L) / L
    return T * H * np.exp(-IntExp) / L

def Permittivity(T):
    D=917
    Freq=1.413e9
    l=299792458/Freq
    e_ice=np.zeros((2,np.size(T)))
    e_ice[0,:] = 3.1884 + 9.1e-4 * (T - 273.0)
    theta = 300.0 / T - 1.0
    alpha = (0.00504 + 0.0062 * theta) * np.exp(-22.1 * theta)
    B1 = 0.0207
    B2 = 1.16e-11
    b = 335
    deltabeta = np.exp(-9.963 + 0.0372 * (T - 273.16))
    betam = (B1 / T) * (np.exp(b / T) / ((np.exp(b / T) - 1) ** 2)) + B2 * (Freq/1e9) ** 2
    beta = betam + deltabeta
    e_ice[1, :] = (alpha /(Freq/1e9) + beta * (Freq/1e9))

    return e_ice

def PlotEmiss(E):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    norm = mpl.colors.Normalize(vmin=0.9, vmax=1.)
    cmap = mpl.cm.spectral
    myplot = ax.pcolormesh(E, cmap=cmap)#, norm=norm)
    cbar = fig.colorbar(myplot, ticks=np.arange(0.9, 1.01, 0.02))
    cbar.set_label('Emissivity', rotation=270)
    cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
    plt.savefig("../../OutputData/img/InvertingEmissDepth/Emissivity_DescentGrad_nDim.png")
    plt.show()

###################################################################################################
#Here load data and compute emissivity
###################################################################################################
Start=time.time()

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
Tb=Tb[0]
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
nc_obsattrs, nc_obsdims, nc_obsvars = ncr.ncdump(Obs)

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4Ts.nc')
H = np.array(GRISLI.variables['H'])
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]
Tz_gr=np.array(Tz_gr)+273.15

Mask=np.zeros(np.shape(H))
Mask[H>=10]=1.0
TbObs = Tb * Mask
Teff=ComputeTeff()
Emissivity = np.ones(np.shape(Teff))
Emissivity[Mask == 1] = (TbObs / Teff)[Mask == 1]
Mask[Emissivity == -32768.0] = 0
Emissivity[Mask == 0] = 0

#PlotEmiss(Mask)
Emissivity=Emissivity*Mask

#PlotEmiss(Teff)
PlotEmiss(Emissivity)

###################################################################################################
#Data output
###################################################################################################

#Create NetCDF file
cols = len(X[0,:])
rows = len(Y[:,0])

outfile = r'../../SourceData/WorkingFiles/Emissivity_FromMatzler_GaussLegendre.nc'
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
Xout[:] = X[0,:]-25000

nc_new.createVariable("Emissivity", 'float64', ('y','x'))
nc_new.variables["Emissivity"][:] = Emissivity[::-1, :]
nc_new.createVariable("Teff", 'float64', ('y','x'))
nc_new.variables["Teff"][:] = Teff[::-1, :]
nc_new.createVariable("Error", 'float64', ('y','x'))
nc_new.variables["Error"][:] = Emissivity[::-1,:]*Teff[::-1, :]-Tb[::-1,:]
crs = nc_new.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
nc_new.close()
