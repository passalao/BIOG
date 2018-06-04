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

def ComputeCorr(Tz_gr,H,Depth,TsMin, TsMax,Subsample):
    print("Depth:", Depth)
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

    #Compute covariance
    m=np.stack(((Emissivity-np.mean(Emissivity))/np.std(Emissivity), (Teff-np.mean(Teff))/np.std(Teff)), axis=0)
    Cov=np.cov(m)
    print("Covariance: ", Cov[0,1])

    '''#if Depth==700:
    cmap = plt.cm.get_cmap("coolwarm")
    norm = mpl.colors.Normalize(vmin=0.95, vmax=1.05)
    plt.plot(x_test, regr.predict(x_test[:, np.newaxis]), color='blue', linewidth=1)
    #plt.plot(x_test, regr.predict(x_test[:, np.newaxis]), color='blue')
    #plt.scatter(TbObs, Teff,s=1, c="r")#, cmap=cmap, norm=norm)
    plt.scatter(TbObs, Teff,s=1, c=Emissivity, cmap=cmap, norm=norm)
    plt.scatter(TbObs, Tbmod,s=0.1, c="g", cmap=cmap, norm=norm)
    plt.plot([210,270],[210,270],c='r')
    plt.xlabel("TbObs")
    plt.ylabel("Tbmod")
    plt.text(220,260,d)
    plt.show()'''
    J1=sum((Tbmod-TbObs)**2)
    #J2=sum((Emissivity-np.mean(Emissivity))**2)
    J2=sum((Emissivity[Emissivity>=1.0]-1)**2)
    return Teff, J1, J2, regr.coef_, r**2, 100*np.std(Emissivity), Cov[0,1]

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
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

TsMax=-50
TsMin=-52.5
Subsample=1
LayerThick=10

#Determine which Depth is the best one
Depths=np.arange(50,1000,25)
#Alpha=np.arange(0,-1e6,-1e4)

#for a in Alpha:
J1 = []
J2 = []
Coeffs = []
Scores=[]
StdEmiss=[]
Cov=[]
for d in Depths:
    Data=ComputeCorr(Tz_gr, H, d, TsMin, TsMax, Subsample)
    J1.append(Data[1])
    J2.append(Data[2])
    Coeffs.append(Data[3])
    Scores.append(Data[4])
    StdEmiss.append(Data[5])
    Cov.append(Data[6])
    #print((max(J1)-min(J1))/(max(J2)-min(J2)))
    #Jtot=J1[-1]+a*J2[-1]
    #print(a,d,Jtot)
print(min(J1), max(J1))
J1=np.array(J1)
J2=np.array(J2)
J1=(J1-min(J1))/(max(J1)-min(J1))
J2=(J2-min(J2))/(max(J2)-min(J2))

Jtot=(J1+J2)/2

plt.plot(Depths,J1,c="b", label="J1: Teff-Tb")
plt.plot(Depths,J2,c="r", label="J2: Emiss-1")
#plt.plot(Depths,Jtot,c="k", label='Total')
plt.plot(Depths,Scores,c="purple", label='Explained variance r$^2$')
plt.plot(Depths,Coeffs,c="orange", label='Mean emissivity')
plt.plot(Depths,StdEmiss,c="green", label='StDev Emissivity')
plt.plot(Depths,Cov,c="pink", label='Cov. Emiss/Teff')
#plt.plot([min(Depths), max(Depths)],[1,1],'--',c='r')
plt.xlabel('Depths (m)')
plt.ylabel('Cost functions, normalized')
plt.legend()
plt.grid()
plt.title(str(TsMin)+ "< Ts < "+str(TsMax))
plt.savefig("../../OutputData/img/InvertingEmissDepth/img.png")
plt.show()