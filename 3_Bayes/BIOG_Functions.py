#!/usr/bin/python
# -*- coding: cp1252 -*-

#  Load the training data that were used for the scaling
#  and exports the corresponding Scaler object
def Load_Scaler(TrainData, TargetLabel):
    import pandas as pd
    from sklearn.preprocessing import StandardScaler
    TrainData = pd.read_csv(TrainData, sep=" ")
    TrainData = TrainData.drop([TargetLabel], axis=1)
    TrainData = TrainData[TrainData.H != 1.0]  # Continent nodes, i.e. where H=1.0
    TrainData = TrainData[TrainData.PhiG < 0.2]
    TrainData = TrainData[TrainData.Acc > 0.]
    MyScaler = StandardScaler()
    MyScaler.fit(TrainData)
    return MyScaler

#Attempt to use theano for inputting suitable data in Neural model within the bayesian inference process...
import theano
import theano.tensor as t
@theano.compile.ops.as_op(itypes=[t.dmatrix,t.dscalar],otypes=[t.dmatrix])
def Update_PhiG(VerticalLine,PhiG):
    VerticalLine[:,3]=PhiG
    print type(VerticalLine)
    return VerticalLine

#Function that leads a MH bayesian inference by exploring a Random Variable (RV) space
import math, random
import BIOG_Classes as bc
def Metropolis(In_Data, RanVarIndexes, Obs_Data, FreeRV_sd, Obs_Data_sd, nbsteps, stepsize, NeuralModel):
    RV=In_Data[0,:]
    #Trace=RV #Start point in the parameter space
    Score=1e6
    Start=time.clock()
    Trace=bc.Trace(np.size(stepsize))

    for n in np.arange(nbsteps):

        if n/100 == float(n)/100:
            print "Step #",n

        PreviousScore=Score
        #Choose new random new value for RV
        i=0
        for rv in RV:
            if np.size(Trace.Data)==np.size(RV):#Concerns the starting point only
                RV[i]=RV[i]+random.uniform(-1,1)*stepsize[i]
            else:
                RV[i]=Trace.Data[-1,i]+random.uniform(-1,1)*stepsize[i]
            i=i+1

        #Consistency control on the data value
        if RV[3]<0: #Geothermal flux should be positive
            pass

        #Evaluate the model at this new point
        #Change the input data with the random value, all along the vertical
        for i in RanVarIndexes:
            In_Data[:,i]=RV[i]
        Tz_mod = NeuralModel.predict(In_Data)
        Tb_mod = GetTb(Tz_mod)

        #Computes the cost function, to be discussed...
        Score=abs(Tb_mod-Obs_Data)
        if np.size(Trace.Data) == np.size(RV):
            Prob = math.exp(-((Tb_mod-Obs_Data)**2/(2*Obs_Data_sd**2)+
                              (RV[3]-Trace.Data[3])**2/(2*FreeRV_sd[0,3]**2)+
                              (RV[2]-Trace.Data[2])**2/(2*FreeRV_sd[0,2]**2)))
        else:
            Prob = math.exp(-((Tb_mod-Obs_Data)**2/(2*Obs_Data_sd**2)+
                              (RV[3]-Trace.Data[-1,3])**2/(2*FreeRV_sd[0,3]**2)+
                              (RV[2]-Trace.Data[-1,2])**2/(2*FreeRV_sd[0,2]**2)))

        #Acception/Rejection process
        RandNum = random.random()
        if Score > PreviousScore :
            if Prob<RandNum:
                if np.size(Trace.Data) == np.size(RV):  #Concerns the starting point only
                    RV=Trace.Data[:]
                else:
                    RV=Trace.Data[-1,:]
            else:
                Trace.Data=np.vstack((Trace.Data,RV))
                Trace.Cost = np.append(Trace.Cost, Prob)
        else:
            Trace.Data=np.vstack((Trace.Data,RV))
            Trace.Cost = np.append(Trace.Cost, Prob)

    Stop = time.clock()
    print "Inference computing time: ", Stop-Start
    #Trace=np.array(Tutu.Trace)
    return Trace


# SMRT model for Tb computation
import sys
sys.path.insert(0, "/home/passalao/smrt")
from smrt import make_snowpack, make_model, sensor
import numpy as np
import time
def GetTb(Tz, H):
    # type: (object, object) -> object
    Start=time.clock()
    # prepare inputs
    l = 20  # number of layers
    if H==0:
        H=1.0#Avoid dummy results

    thickness = np.array([H/l]* l)

    #Prepare temperature field and interpolate on the layer altitude
    temperature = np.transpose(Tz)+273.15
    zini=np.linspace(0,H,21)
    zfin=np.linspace(0,H,l)
    interp_temperature=np.interp(zfin, zini, temperature[::-1])[::-1]

    density=np.zeros(l)
    i=0
    for z in zfin:
       density[i]=1000.*(0.917-0.593*math.exp(-0.01859*z))#From Macelloni et al, 2016
       i=i+1

    p_ex = 1./(917-density+1e-4)
    p_ex[0] = 1e-4

    # create the snowpack
    snowpack = make_snowpack(thickness=thickness,
                             microstructure_model="exponential",
                             density=density,
                             temperature=interp_temperature,
                             corr_length=p_ex)
    # create the snowpack
    m = make_model("iba", "dort")
    # create the sensor
    radiometer = sensor.passive(1.4e9, 52.5)
    # run the model
    res = m.run(radiometer, snowpack)
    # outputs
    return res.TbV()

# Below, functions for K-mean classification
import numpy as np
import random

def cluster_points(X, mu):
    clusters = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x - mu[i[0]])) \
                         for i in enumerate(mu)], key=lambda t: t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters

def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis=0))
    return newmu

def has_converged(mu, oldmu):
    Converge = (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]))
    return Converge

def find_centers(X, K):
    # Initialize to K random centers
    oldmu = random.sample(X, K)
    mu = random.sample(X, K)
    while not has_converged(mu, oldmu):
        oldmu = mu
        # Assign all points in X to clusters
        clusters = cluster_points(X, mu)
        # Reevaluate centers
        mu = reevaluate_centers(oldmu, clusters)
    return (mu, clusters)

def init_board(N):
    X = np.array([(random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)) for i in range(N)])
    return X


import numpy as np


####Fonction de lissage
def Smooth(Lx, Ly, p):
    '''Fonction qui débruite une courbe par une moyenne glissante
    sur 2P+1 points'''
    Lxout = Lx[p: -p]
    Lyout = []
    for index in range(p, len(Ly) - p):
        average = np.mean(Ly[index - p: index + p + 1])
        Lyout.append(average)
    return Lxout, Lyout


#x, y = Smooth(list(range(10)), [2, 4, 6, 8, 6, 5, 4, 5, 6, 8], 1)
#print(x, y)