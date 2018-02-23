#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import pandas as pd
import numpy as np
import netCDF4

# Read in data
print("Import Data")
nc = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLI4KERAS.nc')

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

print("Define train and test datasets")
from sklearn.model_selection import train_test_split
y = T.T #Target label
X = np.concatenate((Zeta.T, H, Acc, Ts, PhiG, LogUh.T, HDiv.T), axis=1) #Random variables

print("Rescale the datasets")
# Standardize the data: rescaling
from sklearn.preprocessing import StandardScaler

# Scale the data with `StandardScaler`
X = StandardScaler().fit_transform(X)

# Split the data up in train and test sets
X_train, X_test, y_train, y_test = train_test_split(pd.DataFrame(X), pd.DataFrame(y), test_size=0.1, random_state=42)

# Create the neural network
print("Now create the neural network")
from keras.models import Sequential
from keras.layers import Dense

model = Sequential()
model.add(Dense(np.shape(X)[1]+1, input_dim=np.shape(X)[1], activation='sigmoid'))
model.add(Dense(np.shape(X)[1]//3*2, activation='sigmoid'))
model.add(Dense(np.shape(X)[1]//3, activation='sigmoid'))
model.add(Dense(21))
model.compile(optimizer='rmsprop', loss='mse', metrics=['mae'])
model.fit(X_train, y_train, epochs=40, verbose=2)
model.save("../../SourceData/WorkingFiles/KERASModels/KERAS_2couches.h5")

y_pred = model.predict(X_test)
mse, mae = model.evaluate(X_test, y_test, verbose=1)
print("score : ", mse ** 0.5, mae)