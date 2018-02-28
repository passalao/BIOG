#!/usr/bin/python
# -*- coding: cp1252 -*-
#  BAYESIAN INFERENCE OF GEOTHERMAL FLUX
#  This script computes posterior distribution of the geothermal flux
#  from GRISLI variable data, emulated by a neural model -> T(z)
#  and SMRT radiative transfer T(z) -> Tb
#  Observation of brightness temperature Tb come from SMOS product

import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG, NC_Resources as ncr
import numpy as np
import netCDF4, time, numbers
from keras.models import load_model
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy

#Import the data
print("Import the model data")
nc = netCDF4.Dataset(BIOG.var.ModelData4Bayes)

H = nc.variables["H"]
LogUh = nc.variables["LogUh"]
PhiG = nc.variables["PhiG"]
T = nc.variables["T"]
Zeta = nc.variables["Zeta"]
Acc = nc.variables["Acc"]
HDiv = nc.variables["HDiv"]
Ts = T[ :, :, 0]
nz=nc.dimensions['z'].size
ny=nc.dimensions['y'].size
nx=nc.dimensions['x'].size
nbXYElts = nx*ny

H_2D=copy.deepcopy(np.array(H))
H = np.reshape(H, nbXYElts)
Acc = np.reshape(Acc, nbXYElts)
PhiG = np.reshape(PhiG, nbXYElts)
Ts = np.reshape(Ts, nbXYElts)

# Modify the dimensionality of the 1D arrays
H.resize(nbXYElts, 1);
Acc.resize(nbXYElts, 1);
PhiG.resize(nbXYElts, 1);
Ts.resize(nbXYElts, 1)
LogUh = np.reshape(LogUh, (nx*ny, nz))
Zeta = np.reshape(Zeta, (nx*ny, nz))
HDiv = np.reshape(HDiv, (nx*ny, nz))
#T = np.reshape(T[:,:, 0:21], (nz,nx*ny))
#Z = Zeta*H

print("Define the dataset and rescale it")
from sklearn.model_selection import train_test_split
X = np.concatenate((Zeta, H, Acc, Ts, PhiG, LogUh, HDiv), axis=1) #Random variables
from sklearn.preprocessing import StandardScaler
Scaler = StandardScaler()
Scaler.fit(X)
X = StandardScaler().fit_transform(X)
ScaledSD=BIOG.var.FreeRV_sd/Scaler.scale_
ScaledSteps=BIOG.var.stepsize/Scaler.scale_
X=np.reshape(X,(ny, nx, 3*nz+4))

print("Import the observation data")
nc_Obs = netCDF4.Dataset(BIOG.var.ObsData)
ny_Obs = nc_Obs.dimensions['cols'].size
nx_Obs = nc_Obs.dimensions['rows'].size
Tb_Obs = nc_Obs.variables['BT_V']
mask = nc_Obs.variables['mask']

print("Load the neural model")
NeuralModel = load_model('../../SourceData/WorkingFiles/KERASModels/KERAS_2couches.h5')

ss=BIOG.var.Subsample

'''# Geographic plot
fig, ax = plt.subplots(nrows=2, ncols=2, sharex='col', sharey='row')
cmap = mpl.cm.seismic
norm = mpl.colors.Normalize(vmin=-10, vmax=10)
myplot=ax[0,0].pcolormesh(T[:,:,0]+273.15-Tb_Obs[0,:,:], cmap=cmap, norm=norm)
#cbar = fig.colorbar(myplot, ticks=np.arange(-10, 10, 1))
#cbar.ax.set_xticklabels(['-10', '0', '10'])  # vertically oriented colorbar
norm = mpl.colors.Normalize(vmin=180, vmax=250)
ax[1,0].pcolormesh(Tb_Obs[0,:,:], cmap=cmap, norm=norm)
ax[1,1].pcolormesh(T[:, :, 0]+273.15-6, cmap=cmap, norm=norm)
ax[0, 0].set_xlim([0, 225.])
ax[0, 0].set_ylim([0, 201.])
ax[1, 1].set_xlim([0, 225.])
ax[1, 1].set_ylim([0, 201.])
#plt.autoscale(True)
#plt.axis('equal')
#plt.savefig("../../OutputData/img/Error_SMOS-sMod_DMRTML.png")
plt.show()
plt.close()'''

for j in np.arange(ny_Obs):
    for i in np.arange(nx_Obs):
        Tb_SMOS = Tb_Obs[0, i, j]
        if mask[i,j]==1.0:#run on the continent only
            if i//ss==i/ss and j//ss==j/ss:
                print("H=",H_2D[i, j])
                #Bayesian inference
                print("Now, Bayesian inference for x=", i, "and y=", j)
                Start=time.clock()
                print("Tb=", Tb_SMOS, "and Ts=", T[i,j,0]+273.15)
                PostData=BIOG.fun.Metropolis(i, j, X, ScaledSteps, BIOG.var.nbsteps, BIOG.var.ColIndex, Tb_SMOS, ScaledSD, BIOG.var.Obs_sd, NeuralModel, Scaler)
                Stop=time.clock()
                print("Time of Bayesian inference: ", Stop-Start, "s")
                UnScaledPostData = Scaler.inverse_transform(PostData.Data)

                #Scatter
                plt.clf()
                #norm = mpl.colors.Normalize(vmin=0, vmax=1)
                plt.scatter(UnScaledPostData[:,24], UnScaledPostData[:,23],  c=np.array(PostData.Cost), s=10., linewidths=0.)
                plt.colorbar()
                plt.grid()
                plt.xlabel("Geothermal flux (mW/m$^2$)")
                plt.ylabel("Surface temperature (K)")
                plt.title("Colorbar: acceptance probability")
                plt.savefig("../../OutputData/img/Bayes/Acceptance_"+str(i)+"_"+str(j)+".png")
                #plt.show()

                '''#Histograms
                plt.clf()
                plt.hist(UnScaledPostData[:,23], normed=True, bins=30)
                plt.xlabel("Ts (K)")
                plt.ylabel("Frequency")
                plt.savefig("../../OutputData/img/Bayes/Hist_Ts_" + str(i) + "_" + str(j) + ".png")

                plt.clf()
                plt.hist(UnScaledPostData[:,24], normed=True, bins=30)
                plt.xlabel("Geothermal Flux (mW/m$^2$)")
                plt.ylabel("Frequency")
                plt.savefig("../../OutputData/img/Bayes/Hist_PhiG_" + str(i) + "_" + str(j) + ".png")
                #Trace
                plt.clf()
                tries=np.arange(len(UnScaledPostData[:,23]))
                plt.plot(tries,np.array(PostData.Cost))
                plt.savefig("../../OutputData/img/Bayes/Trace_" + str(i) + "_" + str(j) + ".png")
                #plt.show()
                plt.close()'''

#TODO... Storing the output data
nc.close()