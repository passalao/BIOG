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

def ComputeJ(Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff = []
    TbObs=[]
    #Tbmod=[]

    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff.append(273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depths[i,j])*0.05*H[i,j]/Depths[i,j])/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depths[i,j])*0.05*H[i,j]/Depths[i,j]))
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
    print("Mean E:", np.mean(Emissivity))
    J3=sum((Emissivity[Emissivity>1]-1)**2)#/np.std(Emissivity)**2

    return Teff, J1, J2, Cov, regr.coef_

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

#Import Slope data
SSlopes=netCDF4.Dataset('../../SourceData/WorkingFiles/SurfaceSlopesOnSMOS.nc')
GradS = SSlopes.variables['Slopes']

TsMax=-50
TsMin=-54
Subsample=1
LayerThick=10

#Determine which Depth is the best one
Depth=np.arange(300,1600,250)
#Emiss=np.arange(0.90,1.01,0.025)
frac=np.arange(0,1.0,0.25)#donne le pourcentage d'éloignement à l'émissivité moyenne
#Si frac=0 =>J2 est minimisé, on minimise J1 ensuite
#Si frac=1 => J1 est minimisé (même nul en fait) et on essaye de minimiser J2

J1 = np.zeros((np.size(Depth), np.size(frac)))
J2 = np.zeros((np.size(Depth), np.size(frac)))
Cov= np.zeros((np.size(Depth), np.size(frac)))
print(np.shape(J1))
u=0
for d in Depth:
    v = 0
    for f in frac:
        print(u,v)
        D=np.ones(np.shape(H))*d
        Data=ComputeJ(Tz_gr, H, D, f, TsMin, TsMax, Subsample)
        J1[u,v]=Data[1]
        J2[u,v]=Data[2]
        Cov[u,v]=Data[3]
        v=v+1
    u = u + 1

J1=np.array(J1)
J2=np.array(J2)

#Compute dJ/dx
def ComputeGrad(J,X, axis):
    dJdx = np.zeros((np.shape(J)[0] - 2, np.shape(J)[1] - 2))
    dx = (max(X) - min(X)) / (np.size(X)-1)
    if axis==0:
        for i in np.arange(0,np.shape(dJdx)[axis],1):
            dJdx[i, :] = (J[i + 1,1:-1] - J[i - 1,1:-1]) / (2*dx)
    if axis==1:
        for i in np.arange(0,np.shape(dJdx)[axis],1):
            dJdx[:,i] = (J[1:-1,i + 1] - J[1:-1,i - 1]) / (2*dx)
            print(J[1:-1,i + 1], J[1:-1,i - 1])
    return dJdx

dJ1dL=ComputeGrad(J1,frac,0)
dJ2dL=ComputeGrad(J2,frac,0)
Lambda1=-dJ1dL/dJ2dL

dJ1de=ComputeGrad(J1,frac,1)
dJ2de=ComputeGrad(J2,frac,1)
Lambda2=-dJ1de/dJ2de
print(Lambda1,Lambda2)

#J1=(J1-min(J1))/(max(J1)-min(J1))
#J2=(J2-min(J2))/(max(J2)-min(J2))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(frac, Depth)
plt.contour(X, Y, J1+Lambda2*J2, 25)#np.arange(0,1e6,2e4))
ax.set_ylabel("Depth")
ax.set_xlabel("Emissivity")
ax.set_zlabel("J1+$\lambda$J2")
plt.grid()
plt.show()