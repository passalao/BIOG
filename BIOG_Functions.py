#!/usr/bin/python
# -*- coding: cp1252 -*-

# SMRT model for Tb computation
import sys
sys.path.insert(0, "/home/passalao/smrt")
from smrt import make_snowpack, make_model, sensor
import numpy as np
import time, math
import sys
sys.path.insert(0, "/home/passalao/dmrtml")
import dmrtml
from pylab import *

def GetTb_DMRTML(Tz, H, NbLayers, Freq, NbStreams, Angle):
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
       radius[i]=1-0.9999*math.exp(-0.01859*z/1e4) #marche bien
       if density[i]>600:
          medium.append('I')
          #radius[i]=1.e-2
       else:
          medium.append('S')
          #radius[i] = 1 - 0.9999 * math.exp(-0.01859 * z / 1e4)  # test
          #radius[i]=1e-4
       i=i+1

    dist = False                  # if True => use RAYLEIGH distribution of particles
    soilp = None # dmrtml.HUTRoughSoilParams(273) # other parameters have their default value
    res = dmrtml.dmrtml(Freq,NbStreams,height,density,radius,temp, tau=dmrtml.NONSTICKY,medium=medium,dist=dist,soilp=soilp)

    CosAngle=math.cos(Angle)
    BestAngle=-9999
    i=0
    for a in res.mhu:
        if abs(a-CosAngle)<abs(BestAngle-CosAngle):
            BestAngle=i
        i=i+1
    #print(BestAngle)
    #print("cosinus: ", np.arccos(res.mhu)/math.pi*180.)
    return res.TbV(BestAngle)



from smrt.permittivity.ice import ice_permittivity_tiuri84
from smrt.permittivity.ice import ice_permittivity_matzler87
def GetTb_SMRT(Tz, H, NbLayers, Freq, Angle, Perm):
    l = NbLayers  # number of layers
    thickness = np.array([H/l]* l)

    #Prepare temperature field and interpolate on the layer altitude
    temperature = np.transpose(Tz)+273.15
    zini=np.linspace(0,H,21)
    zfin=np.linspace(0,H,l)
    interp_temperature=np.interp(zfin, zini, temperature[::-1])[::-1]

    density=np.zeros(l)
    i=0
    for z in zfin:
       density[i]=1000.*(0.916-0.593*math.exp(-0.01859*z))#From Macelloni et al, 2016
       i=i+1

    p_ex = 1./(917-density+100)#100)
    p_ex[0] = 1e-4

    if Perm=="Tiuri":
        ice_permittivity_model = ice_permittivity_tiuri84
    if Perm=="Matzler":
        ice_permittivity_model = ice_permittivity_matzler87

    # create the snowpack
    snowpack = make_snowpack(thickness=thickness,
                             microstructure_model="sticky_hard_spheres",
                             density=density,
                             temperature=interp_temperature,
                             stickiness=0.2, radius=p_ex,
                             corr_length=p_ex,
                             ice_permittivity_model=ice_permittivity_model)

    # create the snowpack
    m = make_model("iba", "dort")
    # create the sensor
    radiometer = sensor.passive(Freq, Angle)
    # run the model
    res = m.run(radiometer, snowpack)
    # outputs
    return res.TbV()

#Function that leads a MH bayesian inference by exploring a Random Variable (RV) space
import math, random, copy
import BIOG_Variables as var
import BIOG_Classes as bc
def Metropolis(xindex, yindex, In_Data, stepsize, nbsteps, RanVarIndexes, Obs_Data, FreeRV_sd, Obs_Data_sd, NeuralModel, Scaler):
    #Add a column to InData to inverse DeltaTb
    nx=np.shape(In_Data)[0]
    ny=np.shape(In_Data)[1]
    #In_Data=np.concatenate(In_Data, np.zeros((nx,ny,1)))

    RV=In_Data[yindex,xindex,:]
    Prior=copy.deepcopy(RV)#In_Data[yindex,xindex,:] #Set the prior value
    Prob=0.
    Start=time.clock()
    Trace=bc.Trace(np.size(stepsize))
    H=Scaler.inverse_transform(RV)[21]#Get the ice thickness

    Bias_ini = 0.0#To correct model bias
    Bias=Bias_ini

    counter=iter(list(np.arange(nbsteps)))
    for n in counter:
        #if n//100 == float(n)/100:
        #   print("Step #",n)

        #Choose new random new value for RV
        i=0
        for rv in RV:
            if np.size(Trace.Data)==np.size(RV):#Concerns the starting point only
                RV[i]=RV[i]+random.uniform(-1,1)*stepsize[i]
            else:
                RV[i]=Trace.Data[-1,i]\
                      +random.uniform(-1,1)*stepsize[i]
            i=i+1

        Bias = Bias + random.uniform(-1,1)*0.0#No bias

        #Check consistency of data
        if Scaler.inverse_transform(RV)[23]>0: #for surface temperature
            RV[23]=(0.0-Scaler.mean_[23])/Scaler.scale_[23]
        if Scaler.inverse_transform(RV)[24]<10: # for geothermal flux
            RV[24]=(10-Scaler.mean_[24])/Scaler.scale_[24]

        #Evaluate the model at this new point
        #Change the input data with the random value
        for i in RanVarIndexes:
            In_Data[yindex, xindex, i]=RV[i]
        ToTest=In_Data[yindex, xindex, :]
        ToTest=np.reshape(ToTest, (1,67))
        Tz_mod = NeuralModel.predict(ToTest)

        if var.RTModel == "DMRT-ML":
            Tb_mod = GetTb_DMRTML(Tz_mod[0], H, var.NbLayers, var.Freq, var.NbStreams, var.Angle)
        if var.RTModel == "SMRT":
            Tb_mod = GetTb_SMRT(Tz_mod[0], H, var.NbLayers, var.Freq, var.Angle, var.Perm)
        Tb_mod = Tb_mod + Bias#Correction of the bias !

        DeltaTb=Tb_mod-Obs_Data

        #Computes the acceptance function, to be discussed...
        PreviousProb = Prob
        Prob = math.exp(-((Tb_mod-Obs_Data)**2/(2*Obs_Data_sd**2)+\
                          (Prior[24]-RV[24])**2/(2*FreeRV_sd[0,24]**2)+\
                          (Prior[23]-RV[23])**2/(2*FreeRV_sd[0,23]**2)))

        #Acception/Rejection process
        RandNum = random.random()
        #print(RandNum, Prob)
        if Prob < PreviousProb:
            if Prob < RandNum:
                if np.size(Trace.Data) == np.size(RV):  #Concerns the starting point only
                    RV=Trace.Data[:]
                    Bias=Bias_ini
                else:
                    RV=Trace.Data[-1,:]
                    Bias=Trace.DeltaTb[-1]
            else:
                Trace.Data=np.vstack((Trace.Data,RV))
                Trace.Cost = np.append(Trace.Cost, Prob)
                Trace.Bias.append(Bias)
                Trace.DeltaTb.append(DeltaTb)
        else:
            Trace.Data=np.vstack((Trace.Data,RV))
            Trace.Cost = np.append(Trace.Cost, Prob)
            Trace.Bias.append(Bias)
            Trace.DeltaTb.append(DeltaTb)
    Stop = time.clock()
    print("Inference computing time: ", Stop-Start)
    return Trace