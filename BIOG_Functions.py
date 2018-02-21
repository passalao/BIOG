#!/usr/bin/python
# -*- coding: cp1252 -*-

# SMRT model for Tb computation
#import sys
#sys.path.insert(0, "/home/passalao/smrt")
#from smrt import make_snowpack, make_model, sensor
import numpy as np
import time, math
import sys
sys.path.insert(0, "/home/passalao/dmrtml")
import dmrtml
from pylab import *

def GetTb_DMRTML(Tz, H, NbLayers, Freq, NbStreams):
    l = NbLayers  # number of layers
    height = np.array([H/l]* l)
    #Prepare temperature field and interpolate on the layer altitude
    temperature = np.transpose(Tz)+273.15
    nbinputlayers=np.size(temperature)
    zini=np.linspace(0,H,nbinputlayers)
    zfin=np.linspace(0,H,l)
    temp=np.interp(zfin, zini, temperature[::-1])[::-1]
    density=np.zeros(l)
    radius=np.zeros(l)
    medium=list()

    i=0
    for z in zfin:
       density[i]=1000.*(0.916-0.593*math.exp(-0.01859*z))#From Macelloni et al, 2016
       #radius[i]=1-0.9999*math.exp(-0.01859*z/1e4) #marche bien
       if density[i]>600:
          medium.append('I')
          radius[i]=1.e-2
       else:
          medium.append('S')
          #radius[i] = 1 - 0.9999 * math.exp(-0.01859 * z / 1e4)  # test
          radius[i]=1e-4
       i=i+1

    dist = False                  # if True => use RAYLEIGH distribution of particles
    soilp = None # dmrtml.HUTRoughSoilParams(273) # other parameters have their default value
    res = dmrtml.dmrtml(Freq,NbStreams,height,density,radius,temp, tau=dmrtml.NONSTICKY,medium=medium,dist=dist,soilp=soilp)
    #print("cosinus: ", res.mhu)
    return res.TbV(0)

#Function that leads a MH bayesian inference by exploring a Random Variable (RV) space
import math, random
import BIOG_Variables as var
import BIOG_Classes as bc
def Metropolis(xindex, yindex, In_Data, stepsize, nbsteps, RanVarIndexes, Obs_Data, FreeRV_sd, Obs_Data_sd, NeuralModel, Scaler):
    RV=In_Data[yindex,xindex,:]
    #RV=In_Data[(var.Xnbcol-1)*y+x,:]
    #Trace=RV #Start point in the parameter space
    Score=1e6
    Start=time.clock()
    Trace=bc.Trace(np.size(stepsize))
    H=Scaler.inverse_transform(RV)[21]#Get the ice thickness

    for n in np.arange(nbsteps):
        if n//100 == float(n)/100:
            print("Step #",n)

        PreviousScore=Score
        #Choose new random new value for RV
        i=0
        for rv in RV:
            if np.size(Trace.Data)==np.size(RV):#Concerns the starting point only
                RV[i]=RV[i]+random.uniform(-1,1)*stepsize[i]
            else:
                RV[i]=Trace.Data[-1,i]\
                      +random.uniform(-1,1)*stepsize[i]
            i=i+1

        #Consistency control on the data value
        if RV[3]<0: #Geothermal flux should be positive
            pass

        #Evaluate the model at this new point
        #Change the input data with the random value, all along the vertical
        for i in RanVarIndexes:
            In_Data[yindex, xindex, i]=RV[i]
        ToTest=In_Data[yindex, xindex, :]
        ToTest=np.reshape(ToTest, (1,67))
        #ToTest.resize(1,67)
        Tz_mod = NeuralModel.predict(ToTest)
        Tb_mod = GetTb_DMRTML(Tz_mod[0], H, var.NbLayers, var.Freq, var.NbStreams)
        #Computes the cost function, to be discussed...
        Score=abs(Tb_mod-Obs_Data)

        if np.size(Trace.Data) == np.size(RV):
            Prob = math.exp(-((Tb_mod-Obs_Data)**2/(2*Obs_Data_sd**2)+
                              (RV[24]-Trace.Data[24])**2/(2*FreeRV_sd[0,24]**2)+
                              (RV[23]-Trace.Data[23])**2/(2*FreeRV_sd[0,23]**2)))
        else:
            Prob = math.exp(-((Tb_mod-Obs_Data)**2/(2*Obs_Data_sd**2)+
                              (RV[24]-Trace.Data[-1,24])**2/(2*FreeRV_sd[0,24]**2)+
                              (RV[23]-Trace.Data[-1,23])**2/(2*FreeRV_sd[0,23]**2)))

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
    print("Inference computing time: ", Stop-Start)
    #Trace=np.array(Tutu.Trace)
    return Trace

'''def GetTb(Tz, H):
    # type: (object, object) -> object
    Start=time.clock()
    # prepare inputs
    l = 10  # number of layers
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
    return res.TbV()'''