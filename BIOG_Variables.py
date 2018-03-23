#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import numpy as np
#Input data and procedures
Data='../../SourceData/WorkingFiles/GRISLI4KERAS.nc'
Xnbcol=387
Ynbcol=374
Znbcol=21
NeuralModelName='../../SourceData/WorkingFiles/KERASModels/KERAS_2couches.h5'
NbAtt=10

#For radiative transfer:
RTModel="DMRT-ML" # SMRT or DMRT-ML
Perm="Matzler"# Matzler or Tiuri for SMRT
NbLayers=1000
Freq=1.4e9 #[Hz] Sensor frequency
NbStreams=16#
Angle=52.5#[deg] View angle of the sensor
Subsample=10# To go faster for tests: 1 point every Subsample points

#For inference process
ModelData4Bayes='../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc'
ObsData='../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc'
stepsize =np.concatenate((np.zeros(21).T, [0.], [0.], [1.], [1.],np.zeros(21).T,np.zeros(21).T))#Step size for each of the variable: Zeta, H, Acc, Ts, PhiG,LogUh, HDiv
FreeRV_sd= np.reshape(np.concatenate((np.zeros(21).T, [0.], [0.], [1.5], [20.], np.zeros(21).T, np.zeros(21).T)),(1,67)) #Standard deviation of the variables, in true physical units
Obs_sd=1.5#Uncertainty on observed Tb
nbsteps = 1000
ColIndex=[23,24]#Indexes of column of the RV that we want to be free