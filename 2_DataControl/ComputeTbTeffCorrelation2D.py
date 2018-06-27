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


def ComputeJ(vars):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff = []
    TbObs=[]
    Depth=vars[0]
    f=vars[1]
    #Tbmod=[]
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
    #J3=sum((Emissivity[Emissivity>1]-1)**2)#/np.std(Emissivity)**2
    #Costs.append([J1,J2,J3])
    #RandomPath.append([Depth, f])

    #Compute explained variance
    xymean=sum((Tbmod-np.mean(Tbmod))*(TbObs-np.mean(TbObs)))/np.size(Tbmod)
    sigmax=(sum((Tbmod-np.mean(Tbmod))**2)/np.size(Tbmod))**0.5
    sigmay=(sum((TbObs-np.mean(TbObs))**2)/np.size(Tbmod))**0.5
    r=xymean/sigmax/sigmay
    #print("explained variance:",r, r**2)

    return J1, J2, (Cov**2)**0.5

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

TsMax=-40
TsMin=-45
Subsample=1
LayerThick=10

'''#Optimisation avec Platypus
RandomPath=[]
Costs=[]
problem = Problem(2, 3)
problem.types[:] = [Real(0, 900), Real(0, 1)]
problem.function = ComputeJ

algorithm = NSGAII(problem)
Start=time.time()
algorithm.run(500)
Stop=time.time()


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#norm = mpl.colors.Normalize(vmin=0.95, vmax=1)
#cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
myplot=ax.scatter([s.objectives[0] for s in algorithm.result],
           [s.objectives[1] for s in algorithm.result],
           [s.objectives[2] for s in algorithm.result])#, c=RandomPath)
#cbar = fig.colorbar(myplot, ticks=np.arange(200, 800, 100))
#cbar.set_label('Depth', rotation=270)
ax.set_xlabel("J1")
ax.set_ylabel("J2")
ax.set_zlabel("J3")
#plt.grid()
plt.show()

plt.plot(np.arange(0,np.size(RandomPath[0])), RandomPath[0])
plt.show()
plt.plot(np.arange(0,np.size(RandomPath[1])), RandomPath[1])
plt.show()'''

#Optimisation manuelle
Depth=np.arange(200,400,25)
frac=np.arange(0.45,0.6,0.05)#donne le pourcentage d'éloignement à l'émissivité moyenne

l1=np.size(Depth)
l2=np.size(frac)
J1 = np.zeros((l1, l2))
J2 = np.zeros((l1, l2))
J3 = np.zeros((l1, l2))
Cov= np.zeros((l1, l2))
print(np.shape(J1))
u=0
for d in Depth:
    v = 0
    for f in frac:
        vars=[d,f]
        print(d,f)
        #D=np.ones(np.shape(H))*d
        Data=ComputeJ(vars)#Tz_gr, H, D, f, TsMin, TsMax, Subsample)
        J1[u,v]=max(-1e316,min(1e316,Data[0]))
        J2[u,v]=max(-1e316,min(1e316,Data[1]))
        J3[u,v]=max(-1e316,min(1e316,Data[2]))
        v=v+1
    u = u + 1

print("J2",J2)
print("J3",J3)

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
    return dJdx

dJ1dL=ComputeGrad(J1,frac,0)
dJ2dL=ComputeGrad(J2,frac,0)
dJ3dL=ComputeGrad(J3,frac,0)

dJ1de=ComputeGrad(J1,frac,1)
dJ2de=ComputeGrad(J2,frac,1)
dJ3de=ComputeGrad(J3,frac,1)
Lambdas=np.zeros(np.shape(dJ1dL))
Mus=np.zeros(np.shape(dJ1dL))

i=0
for u in dJ1dL:
    j=0
    for v in u:
        A=np.matrix([[dJ2de[i,j],dJ3de[i,j]],[dJ2dL[i,j],dJ3dL[i,j]]])
        B=np.matrix([[dJ1de[i,j]],[dJ1dL[i,j]]])
        print(A,np.linalg.inv(A))
        X=np.linalg.solve(A,-B)#np.dot(np.linalg.inv(A),-B)
        Lambdas[i,j]=X[0,0]
        Mus[i, j] = X[1, 0]
        j=j+1
    i=i+1

print("Lambdas", Lambdas)
print("Mus:", Mus)

#Choose the depth so that J3=0 (or J3=min(J3)
'''i,j = np.unravel_index(np.argmin(J3),J3.shape)
m=[k for k, l in enumerate(J2[i,:]) if l == min(J2[i,:])]
Lambda=Lambdas[i-1,m[0]-1]
Mu=Mus[i-1,m[0]-1]'''

Score=Lambdas**2*J2[1:-1,1:-1]**2+Mus**2*J3[1:-1,1:-1]**2
i,m = np.unravel_index(np.argmin(Score),Score.shape)
Lambda=Lambdas[i,m]
Mu=Mus[i,m]

J=J1+Lambda*J2+Mu*J3
u,v = np.unravel_index(np.argmin(J),J.shape)
print("Best solution:", Depth[u], frac[v])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(frac, Depth)
plt.contour(X, Y, J, 25)#J1+3e4* np.arange(0,1e6,2e4)) marche avec 1e4 pour chacun
ax.set_ylabel("Depth")
ax.set_xlabel("Emissivity")
ax.set_zlabel("J")
plt.grid()
plt.show()

