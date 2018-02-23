#!/usr/bin/python
# -*- coding: cp1252 -*-
#  BAYESIAN INFERENCE OF GEOTHERMAL FLUX
#  This script computes posterior distribution of the geothermal flux
#  from GRISLI variable data, emulated by a neural model -> T(z)
#  and SMRT radiative transfer T(z) -> Tb
#  Observation of brightness temperature Tb come from SMOS product

import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import numpy as np
import netCDF4, time
from keras.models import load_model


#Import the data
print("Import the data")
nc = netCDF4.Dataset(BIOG.var.Data)

H = nc.variables["H"]
LogUh = nc.variables["LogUh"]
PhiG = nc.variables["PhiG"]
T = nc.variables["T"]
Zeta = nc.variables["Zeta"]
Acc = nc.variables["Acc"]
HDiv = nc.variables["HDiv"]
Ts = T[0, :, :]
nz=nc.dimensions['z'].size
ny=nc.dimensions['y'].size
nx=nc.dimensions['x'].size
nbXYElts = nx*ny
#print(PhiG[250,250], Ts[250,250])



H = np.reshape(H, nbXYElts)
Acc = np.reshape(Acc, nbXYElts)
PhiG = np.reshape(PhiG, nbXYElts)
Ts = np.reshape(Ts, nbXYElts)

# Modify the dimensionality of the 1D arrays
H.resize(nbXYElts, 1);
Acc.resize(nbXYElts, 1);
PhiG.resize(nbXYElts, 1);
Ts.resize(nbXYElts, 1)
LogUh = np.reshape(LogUh, (nz, nx*ny))
Zeta = np.reshape(Zeta, (nz,nx*ny))
HDiv = np.reshape(HDiv, (nz,nx*ny))
T = np.reshape(T[0:21,:,:], (nz,nx*ny))
Z = Zeta*H.T

print("Define the dataset and rescale it")
from sklearn.model_selection import train_test_split
X = np.concatenate((Zeta.T, H, Acc, Ts, PhiG, LogUh.T, HDiv.T), axis=1) #Random variables
from sklearn.preprocessing import StandardScaler
Scaler = StandardScaler()
Scaler.fit(X)
X = StandardScaler().fit_transform(X)
ScaledSD=BIOG.var.FreeRV_sd/Scaler.scale_
ScaledSteps=BIOG.var.stepsize/Scaler.scale_
X=np.reshape(X,(ny, nx, 3*nz+4))

print("Load the neural model")
NeuralModel = load_model('../../SourceData/WorkingFiles/KERASModels/KERAS_2couches.h5')

#for j in np.arange(ny):
#    for i in np.arange(nx):
i=250; j=250 #Horizontal coordinates where the metropolis will be run
#Bayesian inference
print("Now, Bayesian inference")
Start=time.clock()
Tb_SMOS=220.#Dummy, to be replaced by real SMOS data
PostData=BIOG.fun.Metropolis(i, j, X, ScaledSteps, BIOG.var.nbsteps, BIOG.var.ColIndex, Tb_SMOS, ScaledSD, BIOG.var.Obs_sd, NeuralModel, Scaler)
Stop=time.clock()
print("Time of Bayesian inference: ", Stop-Start, "s")
UnScaledPostData = Scaler.inverse_transform(PostData.Data)

#Plots
import matplotlib.pyplot as plt
import numpy as np
#Scatter
plt.scatter(UnScaledPostData[:,24], UnScaledPostData[:,23], c=np.array(PostData.Cost), s=10., linewidths=0.)
plt.colorbar()
plt.grid()
plt.xlabel("Geothermal flux (mW/m$^2$)")
plt.ylabel("Surface temperature (K)")
plt.title("Colorbar: acceptance probability")
plt.show()

#Histogram
plt.hist(UnScaledPostData[:,23], normed=True, bins=30)
plt.show()
plt.hist(UnScaledPostData[:,24], normed=True, bins=30)
plt.show()
#Trace
tries=np.arange(len(UnScaledPostData[:,23]))
plt.plot(tries,np.array(PostData.Cost))
plt.show()

'''import scipy.interpolate
#Contours of acceptance probability
# Set up a regular grid of interpolation points
x=UnScaledPostData[:,3]
y=UnScaledPostData[:,2]
z=PostData.Cost
xi, yi = np.linspace(0.03, 0.07, 100), np.linspace(-70.,-40., 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
print 'x', x, 'y', y, 'z', z
rbf = scipy.interpolate.Rbf(x, y, z, function='quintic')
zi = rbf(xi, yi)
#Draw the contours
levels = np.linspace(0.8,1.,5)
plt.contour(xi, yi, zi, levels, colors='w')
plt.imshow(zi, vmin=0, vmax=1.0, origin='lower',
           extent=[0.03, 0.07, -70., -40.], aspect=1e-3)

#plt.show()'''
plt.close()
#TODO... Storing the output data