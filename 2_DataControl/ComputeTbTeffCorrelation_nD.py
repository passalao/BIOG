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

def ComputeLagComponents(Depth, E):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    '''Teff = []
    TbObs=[]
    Emissivity=[]'''
    Teff=np.zeros(np.shape(E))
    TbObs=np.zeros(np.shape(E))
    Mask=np.zeros(np.shape(E))
    #Emissivity=np.zeros(np.shape(E))

    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)
                #Teff.append(273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth))
                Mask[i,j]=1
                #TbObs.append(Tb[i,j])
                #Emissivity.append(E[i,j])
    TbObs=Tb*Mask
    #Teff=np.array(Teff)
    #TbObs=np.array(TbObs)
    #Emissivity=np.array(Emissivity)
    #Emissivity = Emissivity[TbObs < 1e45]
    #TbObs = TbObs[TbObs < 1e45]

    Tbmod=Emissivity*Teff

    Teff1D=np.reshape(Teff, (1,np.size(Teff)))[0,:]
    Tbmod1D=np.reshape(Tbmod, (1,np.size(Tbmod)))[0,:]
    TbObs1D=np.reshape(TbObs, (1,np.size(TbObs)))[0,:]
    Emiss1D=np.reshape(Emissivity, (1,np.size(Emissivity)))[0,:]

    Test=(Tbmod1D[Tbmod1D!=0]-TbObs1D[Tbmod1D!=0])**2

    J1=np.sum((Tbmod1D[Tbmod1D!=0]-TbObs1D[Tbmod1D!=0])**2)/np.std(TbObs1D[Tbmod1D!=0])**2
    J2=np.sum((Emiss1D[Emiss1D!=0]-np.mean(Emiss1D[Emiss1D!=0]))**2)

    #Compute normalized covariance = correlation
    m=np.stack(((Emiss1D-np.mean(Emiss1D))/np.std(Emiss1D), (Teff1D-np.mean(Teff1D))/np.std(Teff1D)), axis=0)

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

TsMax=-40
TsMin=-45
Subsample=1

#Optimisation manuelle
Depth=np.arange(200,400,25)
frac=np.arange(0.45,0.6,0.05)#donne le pourcentage d'éloignement à l'émissivité moyenne

DeltaLag=1e10
Lag=1e10
Depth=387
Emissivity=np.random.normal(0.97, 0.005, np.size(Ts))
Emissivity=np.reshape(Emissivity,(np.shape(Ts)))
NewEmiss=Emissivity
Lambda=-1e3
Mu=1e6
dL=25
dE=0.005
stepLambda=1e2
stepMu=1e5
stepE=1e-6
stepL=10

#Now gradient descent
while abs(DeltaLag)>1:
    print(Depth, Emissivity[150,150], Lambda, Mu)
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
    dLagdL=dJ1dL+Lambda*dJ2dL+Mu*dJ3dL
    dLdE=2*Teff*(Tbmod-TbObs)+Mu/N*(Teff-np.mean(Teff)*Mask)+(2*Lambda+Mu/N*Tbmod)*(Emissivity-np.mean(Emiss1D)*Mask)

    #TODO : calculer dJde comme pour dJdL, histoire de comparer...
    #print(dJ1)
    #Compute the new free RVs
    #Mu=Mu-J3*stepMu
    #Lambda=Lambda-J2*stepLambda
    Emissivity=Emissivity-dLdE*stepE
    #Depth=Depth-dLagdL*stepL

    Lag=J1+Lambda*J2+Mu*J3
    DeltaLag=Lag-OldLag
    print("Lagrangien:", Lag, J1, J2, J3)

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

Score=Lambdas**2*J2[1:-1,1:-1]**2+Mus**2*J3[1:-1,1:-1]**2
i,m = np.unravel_index(np.argmin(Score),Score.shape)
Lambda=Lambdas[i,m]
Mu=Mus[i,m]

J=J1+Lambda*J2+Mu*J3
u,v = np.unravel_index(np.argmin(J),J.shape)
print("Best solution:", Depth[u], frac[v])

'''fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(frac, Depth)
plt.contour(X, Y, J, 25)#J1+3e4* np.arange(0,1e6,2e4)) marche avec 1e4 pour chacun
ax.set_ylabel("Depth")
ax.set_xlabel("Emissivity")
ax.set_zlabel("J")
plt.grid()
plt.show()'''