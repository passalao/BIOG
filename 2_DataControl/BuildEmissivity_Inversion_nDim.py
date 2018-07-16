from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time, csv
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import NC_Resources as ncr

###################################################################################################
#Functions
###################################################################################################

#Computes the effective temperature
def ComputeTeff(L):
    Teff=np.zeros(np.shape(H))
    Mask=np.zeros(np.shape(H))
    for i in np.arange(0,np.shape(H)[0], 1):
        for j in np.arange(0,np.shape(H)[1], 1):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>1:
                z=(1-Zeta[i,j])*H[i,j]
                dz=0.05*H[i,j]
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-z/L)*dz/L)/sum(np.exp(-z/L)*dz/L)
                Mask[i,j]=1
    return Teff, Mask

#Computes an initial semi-random emissivity field
def InitEmissivity(L, frac):
    Teff=ComputeTeff(L)[0]
    Mask=ComputeTeff(L)[1]
    TbObs=Tb*Mask
    Emissivity=TbObs/Teff
    Rand=np.random.normal(0, 0.001, np.size(Tb))
    Rand=np.reshape(Rand,(np.shape(Tb)))
    Ebarre=np.mean(Emissivity[Emissivity!=0])
    Emissivity=frac*Ebarre+(1-frac)*Emissivity+Rand
    return Emissivity, Mask

#Computes the cost functions
def ComputeLagComponents(L, E):
    Teff=ComputeTeff(L)[0]
    Mask=ComputeTeff(L)[1]
    TbObs=Tb*Mask
    Tbmod=E*Teff

    Teff1D=np.reshape(Teff, (1,np.size(Teff)))[0,:]
    Tbmod1D=np.reshape(Tbmod, (1,np.size(Tbmod)))[0,:]
    TbObs1D=np.reshape(TbObs, (1,np.size(TbObs)))[0,:]
    Emiss1D=np.reshape(E, (1,np.size(Emissivity)))[0,:]

    J1=(np.sum((Tbmod1D[Tbmod1D!=0]-TbObs1D[Tbmod1D!=0])**2)/np.size(Tbmod1D[Tbmod1D!=0]))
    J2=np.sum((Emiss1D[Emiss1D!=0]-1)**2)/np.size(Tbmod1D[Tbmod1D!=0])

    #Compute normalized covariance = correlation
    m=np.stack((Emiss1D, Teff1D), axis=0)
    Sign=np.sign(np.cov(m)[0,1])
    N=np.size(Mask[Mask==1])
    sigmaE=np.std(Emiss1D[Emiss1D!=0])
    sigmaTe=np.std(Teff1D[Emiss1D!=0])
    J3=abs(np.cov(m)[0,1])/sigmaE/sigmaTe

    #Compute gradients along E
    dLagdE=2*Teff*(Tbmod-TbObs)+Sign*Mu/sigmaTe/sigmaE/N*(Teff-np.mean(Teff1D)*Mask)+(2*Lambda/N+Sign*Mu/sigmaTe/sigmaE/N*Tbmod)*(Emissivity-np.mean(Emiss1D)*Mask)

    print("J1", J1, "J3", J3)
    return J1, J2, J3, dLagdE, Teff

def GradientDescent(Depth, Emissivity, DeltaLag, MinDeltaLag, Lag, J1,J3, Lambda, Mu, dE, dL, stepE, stepL):
    while abs(DeltaLag)>MinDeltaLag and DeltaLag<0:
        OldLag=Lag

        #1 Compute the basic components we need for the lagrangian
        Attempt1 = ComputeLagComponents(Depth, Emissivity)
        J1=Attempt1[0]
        J2=Attempt1[1]
        J3=Attempt1[2]
        dLdE=Attempt1[3]
        Teff=Attempt1[4]

        #2 Compute numerically the gradients along L
        Attempt2=ComputeLagComponents(Depth+dL, Emissivity)
        dJ1dL=(Attempt2[0]-J1)/dL
        dJ2dL=(Attempt2[1]-J2)/dL
        dJ3dL = (Attempt2[2] - J3) / dL
        dLagdL=dJ1dL+Lambda*dJ2dL+Mu*dJ3dL

        #Compute the new free RVs
        Emissivity=Emissivity-dLdE*stepE
        Depth=max(1,Depth-dLagdL*stepL)

        #Update Lagrangian
        Lag=J1+Lambda*J2+Mu*J3
        DeltaLag=Lag-OldLag
        print("DeltaLag",DeltaLag)

    return Emissivity, Teff, Depth, J1, J3

###################################################################################################
#Here solve the optimization problem
###################################################################################################
Start=time.time()

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
Tb=Tb[0]
nc_obsattrs, nc_obsdims, nc_obsvars = ncr.ncdump(Obs)

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4Ts.nc')
H = np.array(GRISLI.variables['H'])
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

#Parameters for data analyse
DeltaT=5
SliceT=np.arange(-60,0,DeltaT)
print("Temperatures slices : ", SliceT)

#For outputs
FinalEmissivity=np.zeros(np.shape(Ts))
FinalTeff=np.zeros(np.shape(Ts))
J1s=[]
J3s=[]
Mus=[]
Depths=[]

#Work slice by slice
for t in SliceT:
    TsMax = t + DeltaT
    TsMin = t

    print(' ')
    print("Slice between ", TsMin, " and ", TsMax)

    MinDeltaLag=1e-2
    MinRelChange=1e-3
    DeltaLag=-1e10
    Lag=1e10
    J1=1e10
    J3=1e10
    Depth=-t*10

    Init=InitEmissivity(Depth, 0.9)
    Emissivity=Init[0]
    Mask=Init[1]
    NewEmiss=Emissivity

    #Initiate gradient descent parameters
    Lambda=0#1e7
    Mu=10
    dL=25
    dE=0.005
    if t<-40:
        stepE=1e-5
    else:
        stepE = 1e-6
    stepL=100000

    #Now gradient descent to optimize L=J1+Lambda*J2+Mu*J3
    Solve=GradientDescent(Depth, Emissivity, DeltaLag, MinDeltaLag, Lag, J1,J3, Lambda, Mu, dE, dL, stepE, stepL)
    Emissivity=Solve[0]
    Teff=Solve[1]
    L=Solve[2]
    J1=Solve[3]
    J3=Solve[4]

    Depths.append(L)
    J1s.append(J1)
    J3s.append(J3)

    FinalEmissivity[Mask==1]=FinalEmissivity[Mask==1]+Emissivity[Mask==1]
    FinalTeff[Mask==1]=FinalTeff[Mask==1]+Teff[Mask==1]

print("J1:", J1s)
print("J3:", J3s)
print("L:", Depths)

Stop=Start=time.time()
print("Elapsed time", Stop-Start, "s")

#Create new NetCDF file for Emissivity data
nc_new = Dataset('../../SourceData/WorkingFiles/Emissivity_FromGradientDescent_nDim.nc', 'w', format='NETCDF4')
nc_new.description = "Emissivity output from inversion of brightness temperature signal"

#Create NetCDF dimensions
for dim in nc_obsdims:
    dim_name=dim
    if dim=="rows":
        dim_name="x"
    if dim=="cols":
        dim_name="y"
    nc_new.createDimension(dim_name, Obs.dimensions[dim].size)
nc_new.createVariable("Emissivity", 'float64', ('x','y'))
nc_new.variables["Emissivity"][:] = FinalEmissivity[:, :]
nc_new.createVariable("Teff", 'float64', ('x','y'))
nc_new.variables["Teff"][:] = FinalTeff[:, :]
nc_new.close()

fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0.9, vmax=1)
cmap = mpl.cm.spectral
myplot = ax.pcolormesh(FinalEmissivity, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0.90, 1.01, 0.01))
cbar.set_label('Emissivity', rotation=270)
cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
plt.savefig("../../OutputData/img/InvertingEmissDepth/Emissivity_DescentGrad_nDim.png")
plt.show()