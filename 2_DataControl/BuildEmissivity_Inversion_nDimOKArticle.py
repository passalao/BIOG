from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time, csv
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import NC_Resources as ncr
import scipy.optimize as opt
#np.set_printoptions(threshold=np.nan)

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
    #Emissivity=TbObs/Teff
    Emissivity=np.ones(np.shape(Teff))
    Emissivity[Mask==1]=(TbObs/Teff)[Mask==1]
    '''i=0
    for e in np.reshape(Emissivity, (1,np.size(Emissivity)))[0,:]:
        if e==-32768:
            print(np.reshape(TbObs, (1,np.size(TbObs)))[0,i])
            print(np.reshape(Mask, (1,np.size(Mask)))[0,i])
        i=i+1'''

    Mask[Emissivity == -32768.0] = 0
    Emissivity[Mask==0]=0

    Rand=np.random.normal(0, 0.00, np.size(Tb))
    Rand=np.reshape(Rand,(np.shape(Tb)))
    Ebarre=np.mean(Emissivity[Emissivity!=0])
    Emissivity=frac*Ebarre+(1-frac)*Emissivity+Rand
    return Emissivity, Mask

#Computes the cost functions
def ComputeLagComponents(L, E, Mask):
    Teff=ComputeTeff(L)[0]
    #Mask=ComputeTeff(L)[1]
    TbObs=Tb*Mask
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
    Sign=np.sign(np.cov(m)[0,1])

    sigmaE=np.std(Emiss1D[Mask1D==1])
    sigmaTe=np.std(Teff1D[Mask1D==1])
    #J3=abs(np.cov(m)[0,1])/sigmaE/sigmaTe
    J3 =((np.cov(m)[0, 1])/sigmaE / sigmaTe)**2

    #Compute gradients along E
    #dLagdE=(2*Teff/N*(Tbmod-TbObs)+2*Lambda/N*(Emissivity-np.mean(Emiss1D)*Mask)+Sign*Mu/sigmaTe/sigmaE/N*(Teff-np.mean(Teff1D)*Mask))*Mask
    dLagdE=(2*Teff/N*(Tbmod-TbObs)+2*Lambda/N*(E-np.mean(Emiss1D[Mask1D==1]))+Mu/N/sigmaTe/sigmaE*J3**0.5*(E-np.mean(Emiss1D[Mask1D==1])))*Mask
    print((2*Teff/N*(Tbmod-TbObs)*Mask)[150,150], (Mu/N/sigmaTe/sigmaE*J3**0.5*(E-np.mean(Emiss1D[Mask1D==1])))[150,150])
    print("J1", J1, "J2", J2, "J3", J3)
    return J1, J2, J3, dLagdE, Teff

def GradientDescent(Depth, Emissivity, DeltaLag, MinDeltaLag, Lag, J1,J3, Lambda, Mu, dL, stepE, stepL, Mask):
    while abs(DeltaLag)>MinDeltaLag and DeltaLag<0:
        OldLag=Lag
        print(Depth, np.mean(Emissivity[Mask==1]), Mu)
        #1 Compute the basic components we need for the lagrangian
        Attempt1 = ComputeLagComponents(Depth, Emissivity, Mask)
        J1=Attempt1[0]
        J2=Attempt1[1]
        J3=Attempt1[2]
        dLdE=Attempt1[3]
        Teff=Attempt1[4]
        #Mu=Mu+J3*stepMu

        #2 Compute numerically the gradients along L
        Attempt2=ComputeLagComponents(Depth+dL, Emissivity, Mask)
        dJ1dL=(Attempt2[0]-J1)/dL
        dJ2dL=(Attempt2[1]-J2)/dL
        dJ3dL=(Attempt2[2]-J3)/dL

        dLagdL=dJ1dL+Lambda*dJ2dL+Mu*dJ3dL

        #Compute the new free RVs
        Emissivity=Emissivity-dLdE*stepE
        Emissivity[Emissivity>1]=1
        Depth=max(1,Depth-dLagdL*stepL)

        #Update Lagrangian
        Lag=J1+Lambda*J2+Mu*J3
        DeltaLag=Lag-OldLag
        print("DeltaLag",DeltaLag)

    return Emissivity, Teff, Depth, J1, J2, J3

def PlotEmiss(E):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    norm = mpl.colors.Normalize(vmin=0.9, vmax=1.0)
    cmap = mpl.cm.spectral
    myplot = ax.pcolormesh(E, cmap=cmap, norm=norm)
    cbar = fig.colorbar(myplot, ticks=np.arange(0.90, 1.01, 0.02))
    cbar.set_label('Emissivity', rotation=270)
    cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
    plt.savefig("../../OutputData/img/InvertingEmissDepth/Emissivity_DescentGrad_nDim.png")
    plt.show()

###################################################################################################
#Data import
###################################################################################################

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


###################################################################################################
#Here solve the optimization problem
###################################################################################################
Start=time.time()

#Parameters for data analyse
DeltaT=5
SliceT=np.arange(-60,0,DeltaT)
print("Temperatures slices : ", SliceT)

#For outputs
FinalEmissivity=np.zeros(np.shape(Ts))
FinalTeff=np.zeros(np.shape(Ts))
J1s=[]
J2s=[]
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
    DeltaLag=-1e10
    Lag=1e10
    J1=1e10
    J2=1e10
    J3=1e10
    #Depth=-t*10
    Depth=90/math.exp(0.036*t) #initiate with plausible depth
    #Could test different value of multiplier to see the best value for each slice

    Init=InitEmissivity(Depth, 0)
    Emissivity=Init[0]
    Mask=Init[1]

    #Initiate gradient descent parameters
    Lambda=0#1e5
    Mu=100
    dL=25
    stepE=5e-6
    stepL=10000
    stepMu=5

    #Now gradient descent to optimize L=J1+Lambda*J2+Mu*J3
    Solve=GradientDescent(Depth, Emissivity, DeltaLag, MinDeltaLag, Lag, J1,J3, Lambda, Mu, dL, stepE, stepL, Mask)
    Emissivity=Solve[0]
    Teff=Solve[1]
    L=Solve[2]
    J1=Solve[3]
    J2=Solve[4]
    J3=Solve[5]

    Depths.append(L)
    J1s.append(J1)
    J2s.append(J2)
    J3s.append(J3)

    FinalEmissivity[Mask==1]=FinalEmissivity[Mask==1]+Emissivity[Mask==1]
    FinalTeff[Mask==1]=FinalTeff[Mask==1]+Teff[Mask==1]

print("J1:", J1s)
print("J2:", J2s)
print("J3:", J3s)
print("L:", Depths)

Stop=time.time()
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
nc_new.createVariable("Error", 'float64', ('x','y'))
nc_new.variables["Error"][:] = FinalEmissivity[:,:]*FinalTeff[:, :]-Tb[:,:]
nc_new.close()

PlotEmiss(FinalEmissivity)