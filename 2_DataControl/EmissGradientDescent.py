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
sys.path.insert(0, "/home/passalao/dmrtml_MLL")
import BIOG
import dmrtml_bis as dmrtml

# m denotes the number of examples here, not the number of features
def gradientDescent(x, y, theta, alpha, m, numIterations, H, TsMin, TsMax, Subsample):
    for i in range(0, numIterations):
        dJd1=0
        dJd2=0

        #Define parameters to optimize
        e0=theta[0]
        e1=theta[1]

        J = ComputeJ(x, H, [e0,e1], TsMin, TsMax, Subsample)
        #J_e1=ComputeJ(x, H, [e0+alpha[0],e1],TsMin, TsMax, Subsample)
        J_e2=ComputeJ(x, H, [e0,e1++alpha[1]],TsMin, TsMax, Subsample)
        #dJde1 = J_e1 - J
        dJde2 = J_e2 - J

        GradJ=np.array([0,dJde2/alpha[1]])
        theta = theta - alpha*GradJ#np.array([dJdl, dJdlam])#+ GradJ*J #alpha*GradJ
        print(J, theta)
    return theta

def Permittivity(Model,T, D, Freq):
    e_ice=np.zeros((2,np.size(T)))

    if Model=="Tiuri":
        e_ice[0,:]=1+1.7*D/1000+0.7*(D/1000)**2
        e_ice[1,:]=1.5871e6*(0.52*D/1000+0.62*(D/1000)**2)*(1/Freq+1.23e-14*Freq**0.5)
        expT=np.exp(0.036*(T-273.15))
        e_ice[1,:]=e_ice[1,:]*expT

    if Model=="Matzler":
        e_ice[0,:] = 3.1884 + 9.1e-4 * (T - 273.0)
        theta = 300.0 / T - 1.0
        alpha = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta)
        B1 = 0.0207
        B2 = 1.16e-11
        b = 335
        deltabeta = exp(-9.963 + 0.0372 * (T - 273.16))
        betam = (B1 / T) * (exp(b / T) / ((exp(b / T) - 1) ** 2)) + B2 * (Freq/1e9) ** 2
        beta = betam + deltabeta
        e_ice[1, :] = alpha /(Freq/1e9) + beta * (Freq/1e9)

    return e_ice

def ComputeJ(Tz_gr,H,Perm,TsMin, TsMax,Subsample):
    e0=Perm[0]
    e1=Perm[1]
    TbObs=[]
    Tbmod=[]
    nbinputlayers = np.size(Tz_gr[0,0,:])
    l = nbinputlayers  # number of layers
    Angle=BIOG.var.Angle
    soilp = dmrtml.HUTRoughSoilParams(273)

    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                thickness = np.array([H[i,j] / l] * l)
                depthini = np.linspace(0, H[i,j], nbinputlayers)
                depthfin = np.linspace(0, H[i,j], l)
                temp=273.15+np.interp(depthfin, depthini,Tz_gr[i,j,:][::-1])[::-1]
                # Compute the permittivity for the whole profile
                e_ice = Permittivity(BIOG.var.Perm, temp, np.ones(21)*917, BIOG.var.Freq)
                res=dmrtml.dmrtml_bis(BIOG.var.Freq, BIOG.var.NbStreams, thickness, 917.,1e-4, temp, medium="I", soilp=soilp, tbatmodown=0, eps_ice=(e_ice[0, :], e0*np.exp((temp-220)*e1)))
                Tbmod.append(res.TbV(Angle))
                TbObs.append(Tb[i,j])

    TbObs=np.array(TbObs)
    Tbmod=np.array(Tbmod)
    Tbmod=Tbmod[TbObs<1e45]
    TbObs = TbObs[TbObs < 1e45]

    J=sum(((Tbmod-TbObs)/np.std(Tbmod))**2)/np.size(Tbmod)

    print(J)

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
#Explore the (epsilon0, epsilon1) space
#where epsilon=epsilon0+T*epsilon1
alpha = [0.2e-7,1e-3] #Steps
theta = [2e-4, 0.05]#2 variables : depth l and lagrangian multiplier lambda
theta = gradientDescent(Tz_gr, Tb, theta, alpha, n, numIterations, H, TsMin, TsMax, Subsample)