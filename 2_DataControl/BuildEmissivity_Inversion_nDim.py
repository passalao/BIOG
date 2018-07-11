from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time, csv
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import NC_Resources as ncr

def InitEmissivity(Depth, frac):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff=np.zeros(np.shape(H))#+1e10
    TbObs=np.zeros(np.shape(H))
    Mask=np.zeros(np.shape(H))

    for i in np.arange(0,np.shape(H)[0], 1):
        for j in np.arange(0,np.shape(H)[1], 1):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)
                Mask[i,j]=1
    TbObs=Tb*Mask
    Emissivity=TbObs/Teff
    #Emissivity[Teff==0]=0
    #print(Emissivity)
    Rand=np.random.normal(0, 0.001, np.size(Tb))
    Rand=np.reshape(Rand,(np.shape(Tb)))
    Ebarre=np.mean(Emissivity[Emissivity!=0])
    Emissivity=frac*Ebarre+(1-frac)*Emissivity+Rand
    return Emissivity, Ebarre

def ComputeLagComponents(Depth, E):
    Teff=np.zeros(np.shape(E))
    Mask=np.zeros(np.shape(E))

    for i in np.arange(0,np.shape(H)[0], 1):
        for j in np.arange(0,np.shape(H)[1], 1):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>1:
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)
                Mask[i,j]=1
    TbObs=Tb*Mask
    Tbmod=E*Teff

    Teff1D=np.reshape(Teff, (1,np.size(Teff)))[0,:]
    Tbmod1D=np.reshape(Tbmod, (1,np.size(Tbmod)))[0,:]
    TbObs1D=np.reshape(TbObs, (1,np.size(TbObs)))[0,:]
    Emiss1D=np.reshape(E, (1,np.size(Emissivity)))[0,:]

    J1=(np.sum((Tbmod1D[Tbmod1D!=0]-TbObs1D[Tbmod1D!=0])**2)/np.size(Tbmod1D))#np.std(TbObs1D[Tbmod1D!=0])**2
    J2=np.sum((Emiss1D[Emiss1D!=0]-np.mean(Emiss1D[Emiss1D!=0]))**2)#/np.std(Emiss1D[Emiss1D!=0])**2

    #Compute normalized covariance = correlation
    m=np.stack((Emiss1D, Teff1D), axis=0)
    signJ3=np.sign(np.cov(m)[0,1])
    J3=(abs(np.cov(m)[0,1])/np.std(Emiss1D[Emiss1D!=0])/np.std(Teff1D[Emiss1D!=0]))
    print("J1", J1, "J3", J3)
    return J1, J2, J3, Teff, Tbmod, TbObs, Teff1D, Tbmod1D, TbObs1D, Emiss1D, Mask, signJ3

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
SliceT=np.arange(-55,-50,DeltaT)
print("Temperatures : ", SliceT)

FinalEmissivity=np.zeros(np.shape(Ts))
#Depths=[]
J1s=[]
J3s=[]
Mus=[]
Ebarres=[]
Depths=np.arange(100,600,50)#
#Work slice by slice
#for t in SliceT:
for d in Depths:
    t=-35
    TsMax = t + DeltaT
    TsMin = t

    print(' ')
    print("Slice between ", TsMin, " and ", TsMax)

    MinDeltaJ1=5e-3
    MinDeltaJ3=5e-3
    MinDeltaLag=2e-2

    DeltaLag=1e10
    DeltaJ1=-1e10
    DeltaJ3=-1e10
    Lag=1e10
    J1=1e10
    J3=1e10
    Depth=d#-t*5

    #Emissivity=np.random.normal(0.97, 0.005, np.size(Ts))
    #Emissivity=np.reshape(Emissivity,(np.shape(Ts)))
    Emissivity=InitEmissivity(Depth, 0.9)[0]
    Ebarre=InitEmissivity(Depth, 0.9)[1]
    #Mask=InitEmissivity(Depth, 0.9)[1]

    NewEmiss=Emissivity
    Lambda=0
    Mu=10
    dL=25
    dE=0.005
    stepLambda=1e2
    stepMu=1e4
    stepE=1e-6
    stepL=100000

    #Now gradient descent to optimize L=J1+Lambda*J2+Mu*J3
    #while abs(DeltaJ3)>MinDeltaJ3 or abs(DeltaJ1)>MinDeltaJ1 <0:
    while abs(DeltaLag)>MinDeltaLag:
        print("Depth: ", Depth, "Emissivity:", np.mean(Emissivity))#, "Lambda: ", Lambda, "Mu: ", Mu)
        OldLag=Lag
        OldJ1=J1
        OldJ3 = J3
        Emiss1D=NewEmiss

        #1 Compute the basic components we need for the lagrangian
        Attempt1 = ComputeLagComponents(Depth, Emissivity)
        J1=Attempt1[0]
        J2=Attempt1[1]
        J3=Attempt1[2]
        Teff=Attempt1[3]
        Tbmod=Attempt1[4]
        TbObs = Attempt1[5]
        Teff1D=Attempt1[6]
        Tbmod1D=Attempt1[7]
        TbObs1D = Attempt1[8]
        Emiss1D = Attempt1[9]
        Mask = Attempt1[10]
        Sign= Attempt1[11]

        N=np.size(Mask[Mask==1])

        # Update Lagrange multipliers
        #Mu = Mu - J3 * stepMu
        Lambda = 0#Lambda - J2 * stepLambda

        #2 Compute the gradients
        Attempt2=ComputeLagComponents(Depth+dL, Emissivity)
        dJ1dL=(Attempt2[0]-J1)/dL
        dJ2dL=(Attempt2[1]-J2)/dL
        dJ3dL = (Attempt2[2] - J3) / dL
        dLagdL=dJ1dL+Lambda*dJ2dL+Mu*dJ3dL

        dLdE=2*Teff*(Tbmod-TbObs)+Sign*Mu/N*(Teff-np.mean(Teff1D)*Mask)+(2*Lambda/N+Sign*Mu/N*Tbmod)*(Emissivity-np.mean(Emiss1D)*Mask)

        dJ1dE=2*Teff*(Tbmod-TbObs)
        dJ2dE=2*Lambda/N*(Emissivity-np.mean(Emiss1D)*Mask)
        dJ3dE=Sign*Mu/N*(Teff-np.mean(Teff1D)*Mask+Sign*(Emissivity-np.mean(Emiss1D)*Mask)*Tbmod)

        #Compute the new free RVs
        Emissivity=Emissivity-dLdE*stepE
        #Emissivity = Emissivity - dJ3dE * stepE
        Depth=max(1,Depth-dLagdL*stepL)

        #New Lagrangian differnce with previous value
        Lag=J1+Lambda*J2+Mu*J3
        DeltaLag=Lag-OldLag
        DeltaJ1=J1-OldJ1
        DeltaJ3=J3-OldJ3
        print("DeltaLag",DeltaLag)

    print("Lambda", Lambda)
    print("Mu:", Mu)
    #Depths.append(Depth)
    Ebarres.append(Ebarre)
    J1s.append(J1)
    J3s.append(J3)

    FinalEmissivity[Mask==1]=FinalEmissivity[Mask==1]+Emissivity[Mask==1]

print("Depths:" ,Depths)
print("J1:", J1s)
print("J3:", J3s)

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

'''c = csv.write(open("Inversion_Values.csv", "wb"))
c.writerow(Depths)
c.writerow(J1s)
c.writerow(J3s)
c.writerow(Mus)

Outputs=np.stack((Depths, J1s, J3s, Mus), axis=1)

print(Outputs)
for o in Outputs[:]:
    print()
    c.writerow(o)
c.close()'''