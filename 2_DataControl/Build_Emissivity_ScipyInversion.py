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
np.set_printoptions(threshold=np.nan)

###################################################################################################
#Functions
###################################################################################################

#Computes the effective temperature
def ComputeTeff(L):
    Teff=np.zeros(np.shape(H))
    #z = (1-Zeta)
    for i in np.arange(0,np.shape(H)[0], 1):
        for j in np.arange(0,np.shape(H)[1], 1):
            if Mask[i,j]==1:
                z = (1 - Zeta[i, j]) * H[i, j]
                dz=0.05*H[i,j]
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-z/L)*dz/L)/sum(np.exp(-z/L)*dz/L)
    return Teff

#Computes mask corresponding to temperature ranges
def ComputeMask(Zone):
    Mask=np.zeros(np.shape(H))
    for i in np.arange(0,np.shape(H)[0], 1):
        for j in np.arange(0,np.shape(H)[1], 1):
            z = (1 - Zeta[i, j]) * H[i, j]
            Tref=np.mean(Tz_gr[i,j,:][z<=500.])
            if Zone == "East":
                if Tref>TsMin and Tref<TsMax and H[i,j]>1 and j>=112:
                    Mask[i,j]=1
            if Zone == "West":
                if Tref>TsMin and Tref<TsMax and H[i,j]>1 and j<112:
                    Mask[i,j]=1
    return Mask

#Computes an initial semi-random emissivity field
def InitEmissivity(L, frac):
    Teff=ComputeTeff(L)
    TbObs=Tb*Mask
    Emissivity=np.ones(np.shape(Teff))
    Emissivity[Mask==1]=(TbObs/Teff)[Mask==1]
    Mask[Emissivity == -32768.0] = 0
    Emissivity[Mask==0]=0

    Rand=np.random.normal(0, 0.001, np.size(Tb))
    Rand=np.reshape(Rand,(np.shape(Tb)))
    Ebarre=np.mean(Emissivity[Emissivity!=0])
    Emissivity=Mask*(frac*Ebarre+(1-frac)*Emissivity+Rand)
    return Emissivity

def ComputeJac(x):
    L=x[0]
    # Compute numerically the gradients along L
    Attempt1 = ComputeLagComponents(x)
    xdL=x
    xdL[0]=L+dL
    Attempt2 = ComputeLagComponents(xdL)
    dJ1dL = (Attempt2[0] - Attempt1[0]) / dL
    dJ2dL = (Attempt2[1] - Attempt1[1]) / dL
    dJ3dL = (Attempt2[2] - Attempt1[2]) / dL
    dLagdL = (dJ1dL + Lambda * dJ2dL + Mu * dJ3dL)

    # The gradients along the Ei are computed analytically
    dLagdE=Attempt2[3]
    dLagdE=np.reshape(dLagdE,(1,np.size(dLagdE)))
    dLagdx = np.concatenate(([dLagdL],  [Attempt1[2]], dLagdE[0]), axis=0)

    return dLagdx

#Computes the cost functions
def ComputeLagComponents(x):
    L=x[0]
    Mu=x[1]
    E=x[2:]
    E=np.reshape(E,np.shape(Tb))
    print("Depth:", L, "Mu:", Mu, "E:", np.mean(E[E!=0]))
    Teff=ComputeTeff(L)
    TbObs=Mask*Tb
    Tbmod=E*Teff
    Teff1D=np.reshape(Teff, (1,np.size(Teff)))[0,:]
    Tbmod1D=np.reshape(Tbmod, (1,np.size(Tbmod)))[0,:]
    TbObs1D=np.reshape(TbObs, (1,np.size(TbObs)))[0,:]
    Emiss1D=np.reshape(E, (1,np.size(Emissivity)))[0,:]
    Mask1D=np.reshape(Mask, (1,np.size(Teff)))[0,:]

    N=np.size(Mask[Mask==1])

    J1=(np.sum((Tbmod1D[Tbmod1D!=0]-TbObs1D[Tbmod1D!=0])**2)/N)#np.size(Tbmod1D[Tbmod1D!=0]))
    J2=np.sum((Emiss1D[Emiss1D!=0]-1)**2)/N#np.size(Tbmod1D[Tbmod1D!=0])

    #Compute normalized covariance = correlation
    m=np.stack((Emiss1D[Mask1D==1], Teff1D[Mask1D==1]), axis=0)
    sigmaE=np.std(Emiss1D[Mask1D==1])
    sigmaTe=np.std(Teff1D[Mask1D==1])
    J3 =((np.cov(m)[0, 1])/sigmaE / sigmaTe)**2

    #Compute gradients along E
    dLagdE=(2*Teff/N*(Tbmod-TbObs)+2*Lambda/N*(E-np.mean(Emiss1D[Mask1D==1]))+2*Mu/N/sigmaTe/sigmaE*J3**0.5*(E-np.mean(Emiss1D[Mask1D==1])))*Mask
    dJ1dE=(2*Teff/N*(Tbmod-TbObs))*Mask
    dJ3dE=(2/N/sigmaTe/sigmaE*J3**0.5*(E-np.mean(Emiss1D[Mask1D==1])))*Mask

    return J1, J2, J3, dLagdE, Teff, dJ1dE, dJ3dE

def ComputeLagrangian(x):
    Solve=ComputeLagComponents(x)
    J1=Solve[0]
    J2=Solve[1]
    J3=Solve[2]
    print("J1", J1, "J3", J3)
    return J1+Mu*J3

def PlotEmiss(E):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    norm = mpl.colors.Normalize(vmin=0.9, vmax=1)
    cmap = mpl.cm.spectral
    myplot = ax.pcolormesh(E, cmap=cmap, norm=norm)
    cbar = fig.colorbar(myplot, ticks=np.arange(0.90, 1.01, 0.02))
    cbar.set_label('Emissivity', rotation=270)
    cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
    plt.savefig("../../OutputData/img/InvertingEmissDepth/Emissivity_DescentGrad_nDim.png")
    plt.show()

###################################################################################################
#Here solve the optimization problem
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

#Parameters for data analyse
#DeltaT=5
#SliceT=np.arange(-60,-15,DeltaT)
#Temperatures=[-60,-50,-45,-40,-35,-30,-25]#for East
#Temperatures=[-40,-35,-30,-25,-20,-15,-10]#for West

#print("Temperatures slices : ", SliceT)

#For outputs
FinalEmissivity=np.zeros(np.shape(Ts))
FinalTeff=np.zeros(np.shape(Ts))

#Work slice by slice
'''for t in SliceT:
    TsMax = t + DeltaT
    TsMin = t'''

Zones=["East"]#, "West"]

for z in Zones:
    if z=="East":
        Temperatures = [-60, -55, -52.5, -50, -45, -40, -35, -30, -25]
    if z=="West":
        Temperatures = [-45, -40, -35, -30, -25, -20, -15, -10]

    i=0
    for t in Temperatures[0:-1]:
        TsMax = Temperatures[i+1]
        TsMin = t
        i=i+1
        print("  ")
        print("Slice between ", TsMin, " and ", TsMax)

        Depth=30/math.exp(0.036*t) #initiate with plausible depth, between T and M
        Mu=100
        Lambda=0
        Mask=ComputeMask(z)
        Emissivity=InitEmissivity(Depth, 0)
        dL=10

        Emissivity=np.reshape(Emissivity,(1,np.size(Emissivity)))
        Bounds=[(0,1)]*np.size(Emissivity)
        Bounds.append((1e-3, 1000))
        Bounds.append((1,1000))
        Bounds.reverse()

        x0=np.concatenate(([Depth], [Mu], Emissivity[0]), axis=0)
        BestValue=opt.fmin_l_bfgs_b(ComputeLagrangian,x0, fprime=ComputeJac, bounds=Bounds, maxiter=200)

        print("Depth:",BestValue[0][0])

        Emissout=np.reshape(BestValue[0][2:], np.shape(Mask))
        FinalEmissivity=FinalEmissivity+Mask*Emissout

#for display
FinalEmissivity[FinalEmissivity==0]=FinalEmissivity[112,100]
FinalEmissivity[FinalEmissivity>1]=1

Stop=time.time()
print("Elapsed time", Stop-Start, "s")

PlotEmiss(FinalEmissivity)

###################################################################################################
#Data output
###################################################################################################

#Create NetCDF file
cols = len(X[0,:])
rows = len(Y[:,0])

outfile = r'../../SourceData/WorkingFiles/Emissivity_FromGradientDescent_Scipy4QGIS.nc'
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
nc_new.variables["Emissivity"][:] = FinalEmissivity[::-1, :]
nc_new.createVariable("Teff", 'float64', ('y','x'))
nc_new.variables["Teff"][:] = FinalTeff[::-1, :]
nc_new.createVariable("Error", 'float64', ('y','x'))
nc_new.variables["Error"][:] = FinalEmissivity[::-1,:]*FinalTeff[::-1, :]-Tb[::-1,:]
crs = nc_new.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
nc_new.close()
