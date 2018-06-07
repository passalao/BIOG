from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
from sklearn import linear_model
import scipy.integrate as integrate
import scipy.special as special
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

# m denotes the number of examples here, not the number of features
def gradientDescent(x, y, theta, alpha, m, numIterations, H, TsMin, TsMax, Subsample):
    for i in range(0, numIterations):
        dJdl=0
        dJdlam=0

        #Define parameters to optimize
        lam=theta[1]
        Depth=theta[0]

        J = ComputeJ(x, H, Depth,TsMin, TsMax, Subsample, lam)
        J_l=ComputeJ(x, H, Depth+alpha[0],TsMin, TsMax, Subsample, lam)
        J_lam=ComputeJ(x, H, Depth,TsMin, TsMax, Subsample, lam+alpha[1])
        dJdl = J_l - J
        dJdlam = J_lam - J
        '''if J_l<J:
            dJdl=J_l-J
        else:
            dJdl = J-J_l
        if J_lam<J:
            dJdlam=J_lam-J
        else:
            dJdlam = J - J_lam'''

        GradJ=np.array([dJdl/alpha[0],dJdlam/alpha[1]])
        theta = theta - alpha*GradJ#np.array([dJdl, dJdlam])#+ GradJ*J #alpha*GradJ
        print(J, theta)
    return theta

def ComputeJ(Tz_gr,H,Depth,TsMin, TsMax,Subsample, lam):
    Teff = []
    TbObs=[]
    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff.append(273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth))
                TbObs.append(Tb[i,j])
    Teff=np.array(Teff)
    TbObs=np.array(TbObs)
    Teff=Teff[TbObs < 1e45]
    TbObs = TbObs[TbObs < 1e45]


    #Compute regression between Teff and TbObs
    regr = linear_model.LinearRegression(fit_intercept = False) #fit_intercept = false : to force the intercept to 0
    regr.fit(Teff[:, np.newaxis], TbObs)
    Tbmod=regr.predict(Teff[:, np.newaxis])
    Score=regr.score(Teff[:, np.newaxis], TbObs)

    '''#Compute regression manually
    xymean=sum((Teff-np.mean(Teff))*(TbObs-np.mean(TbObs)))/np.size(Teff)
    sigmax=(sum((Teff-np.mean(Teff))**2)/np.size(Teff))**0.5
    sigmay=(sum((TbObs-np.mean(TbObs))**2)/np.size(Teff))**0.5
    r=xymean/sigmax/sigmay'''

    Emissivity=np.array(TbObs)/np.array(Teff)

    #Compute covariance
    m=np.stack(((Emissivity-np.mean(Emissivity))/np.std(Emissivity), (Teff-np.mean(Teff))/np.std(Teff)), axis=0)
    Cov=np.cov(m)

    J1=sum(((Tbmod-TbObs)/np.std(Tbmod))**2)/np.size(Tbmod)
    #J2=sum((Emissivity-np.mean(Emissivity))**2)
    #[Emissivity>=1.0]
    if np.size(Emissivity[Emissivity>=1.0])!=0:
        J2=sum(((Emissivity[Emissivity>=1.0]-1)/(np.std(Emissivity))**2))/np.size(Emissivity[Emissivity>=1.0])
    else:
        J2=0

    J=J1+lam*Cov[1,0]
    print("Costs : ", Cov[1,0], lam*Cov[1,0], J1)
    return J#, J1, J2, regr.coef_, r**2, 100*np.std(Emissivity), Cov[0,1]

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']
Tb=Tb[0]
n=np.size(Tb)

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4Ts.nc')
H = np.array(GRISLI.variables['H'])
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

TsMax=-40
TsMin=-45
Subsample=1
LayerThick=10

numIterations= 50
alpha = [10000,1e4] #Steps
theta = [400, 10]#2 variables : depth l and lagrangian multiplier lambda
theta = gradientDescent(Tz_gr, Tb, theta, alpha, n, numIterations, H, TsMin, TsMax, Subsample)