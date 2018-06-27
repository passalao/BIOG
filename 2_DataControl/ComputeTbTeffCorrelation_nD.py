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
import pymc3 as pm
import theano.tensor as tt
import theano

'''#return the value of Teff
def ComputeTeff(Depth):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff = []
    Tb_Obs=[]
    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff.append(273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth))
                Tb_Obs.append(Tb[i,j])
    Teff=np.array(Teff)
    Tb_Obs=np.array(Tb_Obs)
    Teff=Teff[Tb_Obs < 1e45]
    Tb_Obs = Tb_Obs[Tb_Obs < 1e45]
    return Teff, Tb_Obs'''

#return the value of Teff
@theano.compile.ops.as_op(itypes=[t.dscalar, t.dscalar, t.dscalar], otypes=[t.dvector])
def proc_test(Depth):
    #### write an input file
    global X1, X2, i

    params = {};
    params['X1'] = X1;
    params['X2'] = X2;

    params['alpha'] = float(alpha);
    params['beta1'] = float(beta1);
    params['beta2'] = float(beta2);

    ###
    with open("input.txt", 'w') as outfile:
        json.dump(params, outfile)

    #### make a system call
    os.system("python dummy.py input.txt");

    # Get the number of function runs each 100 NOT necessary, just checking
    i += 1
    if i % 100 == 0:
        print
        " number of evaluations {}".format(i)

    #### read the output file and return the value
    mu = np.loadtxt("output.txt");

    return (np.array(mu))


def ComputeTeff(Depth):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff = []
    Tb_Obs=[]
    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff.append(273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth))
                Tb_Obs.append(Tb[i,j])
    Teff=np.array(Teff)
    Tb_Obs=np.array(Tb_Obs)
    Teff=Teff[Tb_Obs < 1e45]
    Tb_Obs = Tb_Obs[Tb_Obs < 1e45]
    return Teff, Tb_Obs

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

#Initiate the process
TsMax=-40
TsMin=-45
Lowerbound=200
Upperbound=300
Subsample=1
LayerThick=10
Depth=(Lowerbound+Upperbound)/2
Teff = ComputeTeff(Depth)[0]
TbObs = ComputeTeff(Depth)[1]
NbDOF=np.size(Teff)

with pm.Model() as model:
    # Unobserved Random Variable, with previous estimations of geothermal flux
    #Depth= pm.distributions.continuous.Uniform('Depth',lower=Lowerbound, upper=Upperbound)
    Depth= pm.Normal('Depth',mu=(Lowerbound+Upperbound)/2,sd=200)
    Emiss= pm.distributions.continuous.Uniform('Emissivity',lower=0.95, upper=1.0,shape=NbDOF)
    Teff=ComputeTeff(Depth)[0]

    #Here define the data
    Tb_mod=Emiss*Teff
    Tb = pm.Normal('Tb', mu=Tb_mod, sd=5, observed=TbObs)
    print("Random variables defined")

    #Sampling process of MCMC. Default sampler is NUTS, but other exists
    #Metropolis, HamiltonianMC...
    step=pm.NUTS()
    trace = pm.sample(100, tune=50, step=step) #njobs keyword allows for mutliple chains in parallel

print(pm.gelman_rubin(trace))

#pm.traceplot(trace);
#pm.plot_posterior(trace);
#plt.show()