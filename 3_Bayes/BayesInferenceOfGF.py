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
import BIOG_Classes as bc
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
H = np.reshape(np.transpose(H), nbXYElts)
Acc = np.reshape(np.transpose(Acc), nbXYElts)
PhiG = np.reshape(np.transpose(PhiG), nbXYElts)
Ts = np.reshape(np.transpose(Ts), nbXYElts)

# Modify the dimensionality of the 1D arrays
H.resize(nbXYElts, 1);
Acc.resize(nbXYElts, 1);
PhiG.resize(nbXYElts, 1);
Ts.resize(nbXYElts, 1)
LogUh = np.reshape(LogUh, (nx*ny, nz))
Zeta = np.reshape(Zeta, (nx*ny, nz))
HDiv = np.reshape(HDiv, (nx*ny, nz))

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
X_Ease2 = nc_Obs.variables['x_ease2']
Y_Ease2 = nc_Obs.variables['y_ease2']

mask = nc_Obs.variables['mask']
print("Load the neural model")
NeuralModel = load_model('../../SourceData/WorkingFiles/KERASModels/KERAS_2couches.h5')

ss=BIOG.var.Subsample
OutputRV=bc.StoreArray(nx_Obs, ny_Obs, 7)
#OutputRV[150,150,3]=81.
#print(np.shape(OutputRV))

'''for j in np.arange(ny_Obs):
    for i in np.arange(nx_Obs):'''

i = 0
for x in X_Ease2[0]:
    j = 0
    for y in Y_Ease2[:,0]:
        if x < 1.5875e6 and x > 1.2125e6 and y > 612500 and y < 987500:#around Dome C
            Tb_SMOS = Tb_Obs[0, j, i]
        #if mask[i,j]==1.0:#run on the continent only
            #if i//ss==i/ss and j//ss==j/ss:
            #if j>160 and j<175 and i<75 and i>45:
            print("H=",H_2D[j, i])
            #Bayesian inference
            print("Bayesian inference for x=", x, "and y=", y)
            Start=time.time()
            print("Tb=", Tb_SMOS, "and Ts=", T[j,i,0]+273.15)
            PostData=BIOG.fun.Metropolis(j, i, X, ScaledSteps, BIOG.var.nbsteps, BIOG.var.ColIndex, Tb_SMOS, ScaledSD, BIOG.var.Obs_sd, NeuralModel, Scaler)
            Stop=time.time()
            print("Time of Bayesian inference: ", Stop-Start, "s")
            UnScaledPostData = Scaler.inverse_transform(PostData.Data)

            indexBest=np.where(PostData.Cost == np.max(PostData.Cost))[0][0]
            OutputRV.Data[i,j,0]= UnScaledPostData[indexBest,23]#np.mean(UnScaledPostData[:np.int(BIOG.var.nbsteps/2),23]) #Ts
            OutputRV.Data[i,j,1]= np.std(UnScaledPostData[:np.int(BIOG.var.nbsteps/2),23]) #sigma Ts
            OutputRV.Data[i,j,2]= UnScaledPostData[indexBest,24] #PhiG
            OutputRV.Data[i,j,3]= np.std(UnScaledPostData[:np.int(BIOG.var.nbsteps/2),24]) #sigma PhiG
            OutputRV.Data[i,j,4]= np.mean(PostData.DeltaTb) #
            OutputRV.Data[i,j,5]= PostData.Cost[indexBest] #Ts
            OutputRV.Data[i,j,6]= PostData.Bias[indexBest] #

            '''#Scatter
            plt.clf()
            #norm = mpl.colors.Normalize(vmin=0, vmax=1)
            plt.scatter(UnScaledPostData[:,24], UnScaledPostData[:,23],  c=np.array(PostData.Cost), s=10., linewidths=0.)
            plt.colorbar()
            plt.grid()
            plt.xlabel("Geothermal flux (mW/m$^2$)")
            plt.ylabel("Surface temperature (K)")
            plt.title("Colorbar: acceptance probability")
            #plt.savefig("../../OutputData/img/Bayes/Acceptance_"+str(i)+"_"+str(j)+".png")
            plt.show()'''

            '''#Histograms
            plt.clf()
            plt.hist(PostData.Bias, normed=True, bins=30)
            plt.xlabel("Bias (K)")
            plt.ylabel("Frequency")
            #plt.savefig("../../OutputData/img/Bayes/Hist_Bias_" + str(i) + "_" + str(j) + ".png")
            plt.show()'''

            '''
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
        j = j + 1
    i = i + 1

# Export of the inferred dataset
w_nc = netCDF4.Dataset('../../OutputData/BayesedDataset_DomeC_FreeBias.nc', 'w', format='NETCDF4')
w_nc.description = "GRISLI data processed by a Bayesian inference "

#for dim in nc_dims:
w_nc.createDimension('x', nc.dimensions['x'].size)
w_nc.createDimension('y', nc.dimensions['y'].size)
w_nc.createVariable('Ts',np.float64, ('x', 'y'))
w_nc.variables['Ts'][:] = OutputRV.Data[:,:,0]
w_nc.createVariable('Ts_std',np.float64, ('x', 'y'))
w_nc.variables['Ts_std'][:] = OutputRV.Data[:,:,1]
w_nc.createVariable('PhiG',np.float64, ('x', 'y'))
w_nc.variables['PhiG'][:] = OutputRV.Data[:,:,2]
w_nc.createVariable('PhiG_std',np.float64, ('x', 'y'))
w_nc.variables['PhiG_std'][:] = OutputRV.Data[:,:,3]
w_nc.createVariable('DeltaTb',np.float64, ('x', 'y'))
w_nc.variables['DeltaTb'][:] = OutputRV.Data[:,:,4]
w_nc.createVariable('Cost',np.float64, ('x', 'y'))
w_nc.variables['Cost'][:] = OutputRV.Data[:,:,5]
w_nc.createVariable('Bias',np.float64, ('x', 'y'))
w_nc.variables['Bias'][:] = OutputRV.Data[:,:,6]

w_nc.close()
nc.close()