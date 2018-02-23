#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import pandas as pd
import netCDF4
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from keras.models import load_model
from sklearn.preprocessing import StandardScaler

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

print("Define train and test datasets")
from sklearn.model_selection import train_test_split
y = T.T #Target label
X = np.concatenate((Zeta.T, H, Acc, Ts, PhiG, LogUh.T, HDiv.T), axis=1) #Random variables

print("Rescale the datasets")
# Standardize the data: rescaling
from sklearn.preprocessing import StandardScaler

# Scale the data with `StandardScaler`
X = StandardScaler().fit_transform(X)

#Call the model and complete the Data with prediction and error
model = load_model('../../SourceData/WorkingFiles/KERASModels/KERAS_2couches.h5')
T_pred=model.predict(X)
Error=T_pred.T-T

# scatterplot
#plt.scatter(Error, T.T, color='DarkGreen', s=0.0001);
#plt.show()

#Profil vertical de T
TempPredProfile=np.reshape(T_pred.T,(nz,ny,nx))
TempProfile=np.reshape(T,(nz,ny,nx))
ZProfile=np.reshape(Z,(nz,ny,nx))
Error=np.reshape(Error,(nz,ny,nx))
Error_zmean= Error.mean(0) #Computes the error along the vertical direction

#Plot Antarctic data
fig, ax = plt.subplots()
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=-3, vmax=3)
#myplot=plt.pcolormesh(Error_zmean, cmap=cmap, norm=norm)
myplot=plt.pcolormesh(Error[0,:,:], cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-3, 4, 1))
cbar.ax.set_xticklabels(['-3', '-2', '-1', '0', '1', '2','3'])  #   vertically oriented colorbar
plt.autoscale(True)
plt.axis('equal')
plt.axis([0, 387, 0, 374])
#plt.savefig("../../OutputData/img/Error_NeuralModel_Surface.png")
plt.show()

'''#Plot transect data
#DomeC : x=275, Pole Sud : y=185
fig, ax = plt.subplots()
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=-3, vmax=3)
myplot=plt.pcolormesh(Error[:,185,:], cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(-3, 4, 1))#[-3, -2,-1,0,1,2, 3])
cbar.ax.set_xticklabels(['-3', '-2', '-1', '0', '1', '2','3'])  #   vertically oriented colorbar
plt.autoscale(True)
plt.axis([0,387,0,21])
plt.gca().invert_yaxis()
plt.savefig("../../OutputData/img/Error_NeuralModel_Transect_PoleSud.png")
plt.show()'''

'''for t in np.arange(0,350,50):
    j=0
    c = ['k', 'b', 'r', 'm', 'g', 'c', 'y']
    plt.clf()
    for i in np.arange(0,350,50):
        plt.plot(TempPredProfile[:,t, i], ZProfile[:,t,i], color=c[j], label="Predicted temp.");
        plt.plot(TempProfile[:,t,i], ZProfile[:,t,i], '--', color=c[j], label="GRISLI temp.");
            #plt.text(-50, i, PhiG[0,250,i])
        j=j+1
    plt.xlabel('Temperature (deg. C)')
    plt.ylabel('Z above bed (m)')
    plt.grid()
    name_fig="../../OutputData/img/Temp_profiles_j"+str(t)+".png"
    plt.savefig(name_fig)'''

c = ['k', 'b', 'r', 'm', 'g', 'c', 'y']

'''j=0
t=200
X=[263,341,288,130]; Y=[284,207,61,138]#Z=0, classe=14
#X=[176,329,265,79]; Y=[210,199,65,194]#Z=0, classe=11
#X=[210,260,288,113];Y=[286,150,79,172]#Z=0 classe=5
#X=[226,293,253,114];Y=[254,139,93,172]#Z=10, classe=14
#X=[266,342,127,176];Y=[281,188,148,214]#Z=10, classe=4
#X=[321,320,134,157];Y=[219,92,135,274]#Z=10, classe=11
#X=[264,324,246,110];Y=[259,173,101,176]#Z=20, classe=12
#X=[156,271,329,129];Y=[271,287,112,136]#Z=20, classe=14
#X=[131,317,148,79];Y=[218,284,135,231]#Z=20, classe=10

for x in X:
    y=Y[j]
#if PhiG[0,t,i]<0.2:#Avoid dummy geo fluxes
    plt.plot(TempPredProfile[:,y, x], ZProfile[:,y,x], color=c[j], label="Predicted temp.");
    plt.plot(TempProfile[:,y,x], ZProfile[:,y,x], '--', color=c[j], label="GRISLI temp.");
    #plt.text(-50, i, PhiG[0,250,i])
    j=j+1

#plt.legend(loc='upper left')
#plt.plot(Data['Temp_pred'].values[1500:1520],Data['Z'].values[1500:1520])
plt.show()'''

print("Elapsed time: ", time.clock(), 's')