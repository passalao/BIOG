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

def ComputeCorrReduced(Tz_gr,H,ReducedDepth,TsMin, TsMax,Subsample):
    Depth=ReducedDepth*H
    Data =ComputeCorr(Tz_gr, H, Depth, TsMin, TsMax, Subsample)
    return Data

def ComputeCorr(Tz_gr,H,Depths,TsMin, TsMax,Subsample):
    Teff = []
    TbObs=[]
    Slopes=[]
    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff.append(273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depths[i,j])*0.05*H[i,j]/Depths[i,j])/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depths[i,j])*0.05*H[i,j]/Depths[i,j]))
                TbObs.append(Tb[i,j])
                Slopes.append(GradS[i,j])
    Teff=np.array(Teff)
    TbObs=np.array(TbObs)
    Slopes=np.array(Slopes)
    Teff=Teff[TbObs < 1e45]
    Slopes=Slopes[TbObs < 1e45]
    TbObs = TbObs[TbObs < 1e45]

    #Compute regression manually
    xymean=sum((Teff-np.mean(Teff))*(TbObs-np.mean(TbObs)))/np.size(Teff)
    sigmax=(sum((Teff-np.mean(Teff))**2)/np.size(Teff))**0.5
    sigmay=(sum((TbObs-np.mean(TbObs))**2)/np.size(Teff))**0.5
    r=xymean/sigmax/sigmay
    #print(xymean, sigmax, sigmay,r)

    #Compute regression between Teff and TbObs
    regr = linear_model.LinearRegression(fit_intercept = False) #fit_intercept = false : to force the intercept to 0
    regr.fit(Teff[:, np.newaxis], TbObs)
    Tbmod=regr.predict(Teff[:, np.newaxis])
    Score=regr.score(Teff[:, np.newaxis], TbObs)

    #Tbmod=regr.coef_*Teff+regr.intercept_
    #print(regr.coef_, regr.intercept_, Score)
    x_test = np.linspace(np.min(Teff), np.max(Teff), 100)

    Emissivity=np.array(TbObs)/np.array(Teff)

    #Compute normalized covariance = correlation
    m=np.stack(((Emissivity-np.mean(Emissivity))/np.std(Emissivity), (Teff-np.mean(Teff))/np.std(Teff)), axis=0)
    m1=np.stack(((Emissivity-np.mean(Emissivity))/np.std(Emissivity), (Slopes-np.mean(Slopes))/np.std(Slopes)), axis=0)
    m2=np.stack(((Teff-np.mean(Teff))/np.std(Teff), (Slopes-np.mean(Slopes))/np.std(Slopes)), axis=0)

    Cov=(np.cov(m))[0,1]
    Cov1=(np.cov(m1))[0,1]
    Cov2=(np.cov(m2))[0,1]

    print("Covariance: ", Cov, Cov1, Cov2)

    #Correction for the non independence of variables
    #Emissivity=Emissivity-rEmissTeff*regr.coef_#Tbmod/Teff
    #print(Emissivity[1000], rEmissTeff)
    #if Depth==700:
    cmap = plt.cm.get_cmap("viridis")
    norm = mpl.colors.Normalize(vmin=220.95, vmax=241.0)
    '''plt.plot(x_test, regr.predict(x_test[:, np.newaxis]), color='blue', linewidth=1)
    #plt.scatter(TbObs, Teff,s=1, c="r")#, cmap=cmap, norm=norm)
    plt.scatter(TbObs, Teff,s=1, c=Emissivity, cmap=cmap, norm=norm)
    plt.scatter(TbObs, Tbmod,s=0.1, c="g", cmap=cmap, norm=norm)
    plt.plot([210,270],[210,270],c='r')
    plt.text(220,260,d)
    plt.xlabel("TbObs")
    plt.ylabel("Tbmod")'''

    plt.scatter(Emissivity,Slopes,s=0.5, c=Teff, cmap=cmap, norm=norm)
    plt.show()

    J1=sum((Tbmod-TbObs)**2)
    J2=np.size(Emissivity[Emissivity>1])/np.size(Emissivity)
    #J2=sum((Emissivity[Emissivity>1]-1)**2)
    #J1=sum(((Tbmod-TbObs)/np.std(Tbmod))**2)/np.size(Tbmod)
    #J2=sum((Emissivity-np.mean(Emissivity))**2)
    #[Emissivity>=1.0]
    '''if np.size(Emissivity[Emissivity>=1.0])!=0:
        J2=sum(((Emissivity[Emissivity>=1.0]-1)/(np.std(Emissivity))**2))#/np.size(Emissivity[Emissivity>=1.0])
    else:
        J2=0'''
    return Teff, J1, J2, regr.coef_, r**2, 100*np.std(Emissivity), Cov

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']
Tb=Tb[0]

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4TsandTb.nc')
H = np.array(GRISLI.variables['H'])
S = np.array(GRISLI.variables['S'])
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

#Import Slope data
SSlopes=netCDF4.Dataset('../../SourceData/WorkingFiles/SurfaceSlopesOnSMOS.nc')
GradS = SSlopes.variables['Slopes']

TsMax=-52.5
TsMin=-55
Subsample=1
LayerThick=10

#Determine which Depth is the best one
Depths=np.arange(300,500,50)
ReducedDepths=np.arange(0.01,1,0.025)
#Alpha=np.arange(0,-1e6,-1e4)

#for a in Alpha:
J1 = []
J2 = []
Coeffs = []
Scores=[]
StdEmiss=[]
Cov=[]

'''for rd in ReducedDepths:
    Data=ComputeCorrReduced(Tz_gr, H, rd, TsMin, TsMax, Subsample)
    J1.append(Data[1])
    J2.append(Data[2])
    Coeffs.append(Data[3])
    Scores.append(Data[4])
    StdEmiss.append(Data[5])
    Cov.append(Data[6])'''

for d in Depths:
    D=np.ones(np.shape(H))*d
    Data=ComputeCorr(Tz_gr, H, D, TsMin, TsMax, Subsample)
    J1.append(Data[1])
    J2.append(Data[2])
    Coeffs.append(Data[3])
    Scores.append(Data[4])
    StdEmiss.append(Data[5])
    Cov.append(Data[6])
    #print((max(J1)-min(J1))/(max(J2)-min(J2)))
    #Jtot=J1[-1]+a*J2[-1]
    #print(a,d,Jtot)
#print(min(J1), max(J1))

J1=np.array(J1)
J2=np.array(J2)
J1=(J1-min(J1))/(max(J1)-min(J1))
J2=(J2-min(J2))/(max(J2)-min(J2))

#Jtot=(J1+J2)/2

plt.plot(Depths,J1,c="b", label="J1: Teff-Tb")
plt.plot(Depths,J2,c="r", label="J2: Emiss-1")
plt.plot(Depths,Scores,c="purple", label='Explained variance r$^2$')
#plt.plot(Depths,Coeffs,c="orange", label='Mean emissivity')
plt.plot(Depths,Cov,c="pink", label='Cov. Emiss/Teff')
'''plt.plot(ReducedDepths,J1,c="b", label="J1: Teff-Tb")
plt.plot(ReducedDepths,J2,c="r", label="J2: Emiss-1")
plt.plot(ReducedDepths,Scores,c="purple", label='Explained variance r$^2$')
#plt.plot(Depths,Coeffs,c="orange", label='Mean emissivity')
plt.plot(ReducedDepths,Cov,c="pink", label='Cov. Emiss/Teff')'''
plt.xlabel('Depths (m)')
plt.ylabel('Cost functions, normalized')
plt.legend()
plt.grid()
plt.title(str(TsMin)+ "< Ts < "+str(TsMax))
plt.savefig("../../OutputData/img/InvertingEmissDepth/img.png")
plt.show()