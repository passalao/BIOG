#!/usr/bin/python
# -*- coding: cp1252 -*-

#Input data and procedures
Data="GRISLI_Data.csv"
Xnbcol=387
Ynbcol=374
Znbcol=21
NeuralModelName='KERAS_on_GRISLI_2couches_HDiv.h5'
NbAtt=10

#For inference process
TargetLabel="Temp" #To be
stepsize =[0., 0., 0.1, 0.1, 0., 0., 0.]#Step size for each of the variable:  H, Uh, Ts, PhiG, Z, Acc, HDiv, in scaled units
FreeRV_sd=[[0., 0., 3., 0.02, 0., 0., 0.]] #Standard deviation of the variables, in true physical units
Obs_sd=1.5#
nbsteps = 100
ColIndex=[2,3]#Indexes of column of the free RV
