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
RTModel="SMRT" # SMRT or DMRT-ML
Perm="Matzler"# Matzler or Tiuri for SMRT
NbLayers=10
Freq=1.4e9 #[Hz] Sensor frequency
NbStreams=64# Number of directions for which radiative equation is solved for
Angle=52.5#[deg] View angle of the sensor
Subsample=50 # to go faster for tests

#For inference process
ModelData4Bayes='../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc'
ObsData='../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc'
stepsize =np.concatenate((np.zeros(21).T, [0.], [0.], [1.], [1.],np.zeros(21).T,np.zeros(21).T))#Step size for each of the variable: Zeta, H, Acc, Ts, PhiG,LogUh, HDiv, in scaled units
FreeRV_sd= np.reshape(np.concatenate((np.zeros(21).T, [0.], [0.], [1.5], [20.], np.zeros(21).T, np.zeros(21).T)),(1,67)) #Standard deviation of the variables, in true physical units
Obs_sd=1.5#
nbsteps = 1000
ColIndex=[23,24]#Indexes of column of the free RV