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
import theano.tensor as tt
import pymc3 as pm
import seaborn as sns

sns.set_context('notebook')
plt.style.use('seaborn-darkgrid')
print('Running on PyMC3 v{}'.format(pm.__version__))

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

def ComputeTb(Epsilon0):#Tz,H,Perm):#,TsMin, TsMax,Subsample):
    Thick=H[120,120]
    Tz=Tz_gr[120,120]
    e0=Epsilon0#Perm[0]
    e1=Epsilon1#Perm[1]
    TbObs=[]
    Tbmod=[]
    nbinputlayers = np.size(Tz_gr[0,0,:])
    l = nbinputlayers  # number of layers
    Angle=BIOG.var.Angle
    soilp = dmrtml.HUTRoughSoilParams(273)

    '''for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:'''
    thickness = np.array([Thick / l] * l)
    depthini = np.linspace(0, Thick, nbinputlayers)
    depthfin = np.linspace(0, Thick, l)
    temp=273.15+np.interp(depthfin, depthini,Tz[::-1])[::-1]
    # Compute the permittivity for the whole profile
    e_ice = Permittivity(BIOG.var.Perm, temp, np.ones(21)*917, BIOG.var.Freq)
    res=dmrtml.dmrtml_bis(BIOG.var.Freq, BIOG.var.NbStreams, thickness, 917.,1e-4, temp, medium="I", soilp=soilp, tbatmodown=0, eps_ice=(e_ice[0, :], e0*np.exp((temp-220)*e1)))

    return res.TbV(Angle)#, J1, J2, regr.coef_, r**2, 100*np.std(Emissivity), Cov[0,1]

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

#Explore the (epsilon0, epsilon1) space
#where epsilon=epsilon0+T*epsilon1
theta = [1.2e-4, 0.04]#2 variables epsilon0 and epsilon1
Epsilon1=0.04
Epsilon0=1.2e-4
print(ComputeTb(Epsilon0))#Tz_gr[120,120], H[120,120], theta))
Angle = BIOG.var.Angle
soilp = dmrtml.HUTRoughSoilParams(273)

with pm.Model() as model:
    # Unobserved Random Variable
    Epsilon0 = pm.Normal('Epsilon0', mu=1.2e-4, sd=0.4e-4)#Unobserved RV
    Tb_obs=230+10*np.random.randn(100) #Dummy values for the moment, One observation
    theta=[Epsilon0,Epsilon1]
    print(Epsilon0.logp({'Epsilon0': 99}))

    #Here define the data
    e_ice = Permittivity(BIOG.var.Perm, Tz_gr[120,120], np.ones(21)*917, BIOG.var.Freq)
    res=pm.Deterministic('titi',dmrtml.dmrtml_bis(BIOG.var.Freq, BIOG.var.NbStreams, H[120,120], 917.,1e-4, Tz_gr[120,120], medium="I", soilp=soilp, tbatmodown=0, eps_ice=(e_ice[0, :], Epsilon0*np.exp((Tz_gr[120,120]-220)*Epsilon1))))
    Tb_mod=res.TbV(Angle)#ComputeTb(Epsilon0)#Tz_gr[120,120], H[120,120], theta)
    #Tb_mod=175+10*np.random.randn(100) #KERASmodel.predict(X)+273.15 #Dummy values for the moment, this should be done from the ice+RT model
    Tb = pm.Normal('Tb', mu=Tb_mod, sd=5, observed=Tb_obs)
    print("Random variables defined")

    #Sampling process of MCMC. Default sampler is NUTS, but other exists
    #Metropolis, HamiltonianMC...
    step=pm.NUTS()
    trace = pm.sample(1000, tune=500, step=step) #njobs keyword allows for mutliple chains in parallel

print(pm.gelman_rubin(trace))

#pm.traceplot(trace);
pm.plot_posterior(trace);
plt.show()

