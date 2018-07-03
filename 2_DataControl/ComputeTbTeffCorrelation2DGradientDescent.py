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

def ComputeLagComponents(Depth, f):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff = []
    TbObs=[]

    print(Depth,f)
    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff.append(273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth))
                TbObs.append(Tb[i,j])
    Teff=np.array(Teff)
    TbObs=np.array(TbObs)
    Teff=Teff[TbObs < 1e45]
    TbObs = TbObs[TbObs < 1e45]

    #Compute mean emissivity
    regr = linear_model.LinearRegression(fit_intercept=False)  # fit_intercept = false : to force the intercept to 0
    regr.fit(Teff[:, np.newaxis], TbObs)

    ebound=np.array(TbObs)/np.array(Teff)
    ebarre=np.ones(np.shape(ebound))*regr.coef_
    Emissivity=ebarre*(1-f)+ebound*f
    Tbmod=Emissivity*Teff

    #Compute normalized covariance = correlation
    m=np.stack(((Emissivity-np.mean(Emissivity))/np.std(Emissivity), (Teff-np.mean(Teff))/np.std(Teff)), axis=0)
    Cov=(np.cov(m))[0,1]

    J1=sum((Tbmod-TbObs)**2)/np.std(TbObs)**2
    J2=sum((Emissivity-np.mean(Emissivity))**2)

    return J1, J2, (Cov**2)**0.5, Teff, Tbmod, TbObs, Emissivity

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

TsMax=-50
TsMin=-55
Subsample=1

DeltaLag=1e10
Lag=1e10
Depth=600
Frac=0.5
#Emissivity=0.98#np.random.normal(0.97, 0.005, np.size(Ts))
#Emissivity=np.reshape(Emissivity,(np.shape(Ts)))
#NewEmiss=Emissivity
Lambda=10000
Mu=100
dL=25
df=0.001
stepLambda=1e3
stepMu=1e2
stepf=1e-4
stepL=1e2

#Now gradient descent
while abs(DeltaLag)>2:
    print("Input :", Depth, Frac, Lambda, Mu)
    OldLag=Lag
    #Emiss1D=NewEmiss
    #1 Compute the basic components we need for the lagrangian
    Attempt1 = ComputeLagComponents(Depth, Frac)
    J1=Attempt1[0] #max(-1e316,min(1e316,Attempt1[0]))
    J2=Attempt1[1] #max(-1e316,min(1e316,Attempt1[1]))
    J3=Attempt1[2] #max(-1e316,min(1e316,Attempt1[2]))
    Teff1D = Attempt1[3]
    Tbmod1D = Attempt1[4]
    TbObs1D = Attempt1[5]
    Emiss1D = Attempt1[6]
    #2 Compute the gradients
    Attempt2=ComputeLagComponents(Depth+dL, Frac)
    Attempt3=ComputeLagComponents(Depth, Frac+df)

    dJ1dL=(Attempt2[0]-J1)/dL
    dJ2dL=(Attempt2[1]-J2)/dL
    dJ3dL = (Attempt2[2] - J3) / dL
    dJ1df=(Attempt3[0]-J1)/df
    dJ2df=(Attempt3[1]-J2)/df
    dJ3df = (Attempt3[2] - J3) / df

    dLagdL=dJ1dL+Lambda*dJ2dL+Mu*dJ3dL
    dLdf = dJ1df + Lambda * dJ2df + Mu * dJ3df
    N=np.size(Teff1D)

    #dLdE=2*Teff*(Tbmod-TbObs)+Mu/N*(Teff-np.mean(Teff)*Mask)+(2*Lambda+Mu/N*Tbmod)*(Emissivity-np.mean(Emiss1D)*Mask)
    #Compute the new free RVs
    Mu=Mu-J3*stepMu
    Lambda = Lambda - J2 * stepLambda

    #Pour aider à trouver les valeurs de Lambda et Mu
    #Mu=-dJ1dL/dJ3dL
    #Lambda=np.mean(0.5*(-2/(Emiss1D-np.mean(Emiss1D)*(Teff1D*(Tbmod1D-TbObs1D)-Mu/N*(Teff1D-np.mean(Teff1D))))-Mu/N*Tbmod1D))
    Frac=Frac-dLdf*stepf
    Depth=Depth-dLagdL*stepL

    Lag=J1+Lambda*J2+Mu*J3
    DeltaLag=Lag-OldLag
    print("Lagrangien:", Lag, J1, J2, J3)

print("Lambda", Lambda)
print("Mu:", Mu)

'''fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(frac, Depth)
plt.contour(X, Y, J, 25)#J1+3e4* np.arange(0,1e6,2e4)) marche avec 1e4 pour chacun
ax.set_ylabel("Depth")
ax.set_xlabel("Emissivity")
ax.set_zlabel("J")
plt.grid()
plt.show()'''