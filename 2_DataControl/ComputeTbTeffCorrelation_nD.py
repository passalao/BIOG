from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.stats import chi2_contingency as chi2
from sklearn import linear_model
import scipy.integrate as integrate
import scipy.special as special
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
from platypus import NSGAII, Problem, Real
#np.set_printoptions(threshold=np.nan)

#Compute dJ/dx
def ComputeGrad(J,dx, axis):
    dJdx = np.zeros((np.shape(J)[0] - 2, np.shape(J)[1] - 2))
    #dx = (max(X) - min(X)) / (np.size(X)-1)
    if axis==0:
        for i in np.arange(0,np.shape(dJdx)[axis],1):
            dJdx[i, :] = (J[i + 1,1:-1] - J[i - 1,1:-1]) / (2*dx)
    if axis==1:
        for i in np.arange(0,np.shape(dJdx)[axis],1):
            dJdx[:,i] = (J[1:-1,i + 1] - J[1:-1,i - 1]) / (2*dx)
    return dJdx


def ComputeEmissivity(Depth):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff=np.zeros(np.shape(H))
    TbObs=np.zeros(np.shape(H))
    Mask=np.zeros(np.shape(H))

    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)
                Mask[i,j]=1
    TbObs=Tb*Mask
    Emissivity=TbObs/Teff
    Ebarre=np.mean(Emissivity)
    Emissivity=(Emissivity+Ebarre)/2
    return Emissivity

def ComputeLagComponents(Depth, E):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff=np.zeros(np.shape(E))
    TbObs=np.zeros(np.shape(E))
    Mask=np.zeros(np.shape(E))

    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)
                Mask[i,j]=1
    TbObs=Tb*Mask
    Tbmod=E*Teff

    Teff1D=np.reshape(Teff, (1,np.size(Teff)))[0,:]
    Tbmod1D=np.reshape(Tbmod, (1,np.size(Tbmod)))[0,:]
    TbObs1D=np.reshape(TbObs, (1,np.size(TbObs)))[0,:]
    Emiss1D=np.reshape(Emissivity, (1,np.size(Emissivity)))[0,:]

    Test=(Tbmod1D[Tbmod1D!=0]-TbObs1D[Tbmod1D!=0])**2

    J1=np.sum((Tbmod1D[Tbmod1D!=0]-TbObs1D[Tbmod1D!=0])**2)/np.std(TbObs1D[Tbmod1D!=0])**2
    J2=np.sum((Emiss1D[Emiss1D!=0]-np.mean(Emiss1D[Emiss1D!=0]))**2)#/np.std(Emiss1D[Emiss1D!=0])**2

    #Compute normalized covariance = correlation
    m=np.stack((Emiss1D, Teff1D), axis=0)
    #/np.std(Emiss1D) et /np.std(Teff1D)
    J3=(((np.cov(m))[0,1])**2)**0.5

    return J1, J2, J3, Teff, Tbmod, TbObs, Teff1D, Tbmod1D, TbObs1D, Emiss1D, Mask

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
S = np.array(GRISLI.variables['S'])
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

TsMax=-30
TsMin=-35
Subsample=1

#Optimisation manuelle
Depth=np.arange(260,400,25)
frac=np.arange(0.45,0.6,0.05)#donne le pourcentage d'éloignement à l'émissivité moyenne

DeltaLag=1e10
Lag=1e10
Depth=260
#Emissivity=ComputeEmissivity(Depth)
#Emissivity=Emissivity+np.reshape(np.random.normal(0, 0.005, np.size(Ts)),(np.shape(Ts)))
Emissivity=np.random.normal(0.97, 0.005, np.size(Ts))
Emissivity=np.reshape(Emissivity,(np.shape(Ts)))
NewEmiss=Emissivity
Lambda=0#-115#285#5e3
Mu=0#870959#-2160021#-1e2
dL=25
dE=0.005
stepLambda=1e2
stepMu=1e1
stepE=1e-6
stepL=10

#Now gradient descent
while abs(DeltaLag)>200:
    print(Depth, Lambda, Mu)
    OldLag=Lag
    Emiss1D=NewEmiss
    #1 Compute the basic components we need for the lagrangian
    Attempt1 = ComputeLagComponents(Depth, Emissivity)
    J1=Attempt1[0] #max(-1e316,min(1e316,Attempt1[0]))
    J2=Attempt1[1] #max(-1e316,min(1e316,Attempt1[1]))
    J3=Attempt1[2] #max(-1e316,min(1e316,Attempt1[2]))
    Teff=Attempt1[3]
    Tbmod=Attempt1[4]
    TbObs = Attempt1[5]
    Teff1D=Attempt1[6]
    Tbmod1D=Attempt1[7]
    TbObs1D = Attempt1[8]
    Emiss1D = Attempt1[9]
    Mask = Attempt1[10]

    N=np.size(Emiss1D)

    #2 Compute the gradients
    Attempt2=ComputeLagComponents(Depth+dL, Emissivity)
    dJ1dL=(Attempt2[0]-J1)/dL
    dJ2dL=(Attempt2[1]-J2)/dL
    dJ3dL = (Attempt2[2] - J3) / dL
    '''Attempt3=ComputeLagComponents(Depth, Emissivity+dE)
    dJ1dE=(Attempt3[0]-J1)/dE
    dJ2dE=(Attempt3[1]-J2)/dE
    dJ3dE = (Attempt3[2] - J3) / dE'''

    #Mu = -dJ1dL / dJ3dL
    #Lambda=np.mean(0.5*(-2/(Emiss1D-np.mean(Emiss1D)*(Teff1D*(Tbmod1D-TbObs1D)-Mu/N*(Teff1D-np.mean(Teff1D))))-Mu/N*Tbmod1D))

    #New Lagrange multipliers
    Mu = Mu - J3 * stepMu
    Lambda=Lambda-J2*stepLambda

    #Lagrangian gradients
    dLagdL=dJ1dL+Lambda*dJ2dL+Mu*dJ3dL
    #dLdE = dJ1dE + Lambda * dJ2dE + Mu * dJ3dE
    dLdE=2*Teff*(Tbmod-TbObs)+Mu/N*(Teff-np.mean(Teff1D)*Mask)+(2*Lambda/N+Mu/N*Tbmod)*(Emissivity-np.mean(Emiss1D)*Mask)
    #dLdE=2*Teff*(Tbmod-TbObs)+2*Lambda*(Emissivity-np.mean(Emiss1D)*Mask)
    dJ1dE=2*Teff*(Tbmod-TbObs)
    dJ2dE=2*Lambda/N*(Emissivity-np.mean(Emiss1D)*Mask)
    dJ3dE=Mu/N*(Teff-np.mean(Teff1D)*Mask+(Emissivity-np.mean(Emiss1D)*Mask)*Tbmod)
    print("dJ1dE:", np.mean(dJ1dE))
    print("dJ2dE:", Lambda*np.mean(dJ2dE))
    print("dJ3dE:", Mu*np.mean(dJ3dE))

    #Compute the new free RVs
    Emissivity=Emissivity-dJ1dE*stepE
    Depth=Depth-dLagdL*stepL

    #New Lagrangiana dn differnce with previous value
    Lag=J1+Lambda*J2+Mu*J3
    DeltaLag=Lag-OldLag

    print("Lagrangien:", Lag, J1, J2, J3)
    print("DeltaLag : ", DeltaLag)
print("Lambda", Lambda)
print("Mu:", Mu)

fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0.9, vmax=1)
cmap = mpl.cm.spectral
myplot = ax.pcolormesh(Emissivity, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0.90, 1.01, 0.01))
cbar.set_label('Emissivity', rotation=270)
cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
plt.show()