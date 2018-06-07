from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import numpy as np
import pyproj
from sklearn import linear_model
import scipy.integrate as integrate
import scipy.special as special
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
from platypus import NSGAII, Problem, Real

def ComputeJ(Depth):#Tz_gr,H,Depth,TsMin, TsMax,Subsample):
    RandomPath.append(Depth)
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
    Cov=np.cov(m)[0,1]**2

    J1=sum((Tbmod-TbObs)**2)
    J2=sum((Emissivity[Emissivity>1]-1)**2)

    Costs.append([J1,J2])
    print(Depth, J1, Cov)
    return [J1, J2, Cov]#, Cov]#, regr.coef_, r**2, 100*np.std(Emissivity), Cov[0,1]

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

TsMax=-30
TsMin=-35
Subsample=1
LayerThick=10

RandomPath=[]
Costs=[]

problem = Problem(1, 3)
problem.types[:] = Real(100, 500)
problem.function = ComputeJ

algorithm = NSGAII(problem)
algorithm.run(1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#norm = mpl.colors.Normalize(vmin=0.95, vmax=1)
#cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])

myplot=ax.scatter([s.objectives[0] for s in algorithm.result],
           [s.objectives[1] for s in algorithm.result],
           [s.objectives[2] for s in algorithm.result], c=RandomPath)

#myplot=plt.scatter([s.objectives[0] for s in algorithm.result], [s.objectives[1] for s in algorithm.result], c=RandomPath)
cbar = fig.colorbar(myplot, ticks=np.arange(200, 800, 100))
cbar.set_label('Depth', rotation=270)
ax.set_xlabel("J1")
ax.set_ylabel("J2")
ax.set_zlabel("Cov")
#plt.grid()
plt.show()
