#!/usr/bin/python
# -*- coding: cp1252 -*-
#  BAYESIAN INFERENCE OF GEOTHERMAL FLUX
#  This script computes posterior distribution of the geothermal flux
#  from GRISLI variable data, emulated by a neural model -> T(z)
#  and SMRT radiative transfer T(z) -> Tb
#  Observation of brightness temperature Tb come from SMOS product

import BIOG_Variables as var
import BIOG_Functions as fun
import numpy as np

#  Load the neural model for emulation, and corresponding scaler
import pandas as pd
print "Load neural model"
from keras.models import load_model
NeuralModel = load_model(var.NeuralModelName)
#The scaler should correspond to the dataset used when building the neural model
Scaler=fun.Load_Scaler(var.Data, var.TargetLabel)

#  Get the variable data of GRISLI, isolate the target, and rescale
print "Get the GRISLI data"
TestData=pd.read_csv(var.Data, sep=" ")
GRISLIData = TestData.drop([var.TargetLabel], axis=1)
ScaledData = Scaler.transform(GRISLIData)
ScaledVariances = Scaler.transform(var.FreeRV_sd)

print "Reshaping the data"
#  We need to work vertical line by vertical line
#  Reshape the data to get back the info on verticality, select vertical line
ScaledDataReshape = np.reshape(ScaledData,(var.Znbcol,var.Ynbcol,var.Xnbcol,np.shape(ScaledData)[1]))#
VerticalLine=ScaledDataReshape[:,250,250,:]#  We take only one vertical line. TODO: Replace the figure by incremental indice

#Bayesian inference
import time
Start=time.clock()
print "Now, Bayesian inference"
Tb_SMOS=225.#Dummy, to be replaced by real SMOS data
PostData=fun.Metropolis(VerticalLine, var.ColIndex, Tb_SMOS, ScaledVariances, var.Obs_sd, var.nbsteps, var.stepsize, NeuralModel)
Stop=time.clock()
print "Time of Bayesian inference: ", Stop-Start, "s"
UnScaledPostData = Scaler.inverse_transform(PostData.Data)

#Plots
import matplotlib.pyplot as plt
import numpy as np
#Scatter
plt.scatter(UnScaledPostData[:,3], UnScaledPostData[:,2], c=np.array(PostData.Cost), s=10., linewidths=0.)
plt.colorbar()
plt.show()

'''#Histogram
plt.hist(UnScaledPostData[:,3], normed=True, bins=30)
plt.show()
plt.hist(UnScaledPostData[:,2], normed=True, bins=30)
plt.show()
#Trace
tries=np.arange(len(UnScaledPostData[:,2]))
plt.plot(tries,np.array(PostData.Cost))
plt.show()'''

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

#TODO... Storing the output data