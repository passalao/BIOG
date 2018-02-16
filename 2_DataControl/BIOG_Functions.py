#!/usr/bin/python
# -*- coding: cp1252 -*-

# SMRT model for Tb computation
import sys
sys.path.insert(0, "/home/passalao/smrt")
from smrt import make_snowpack, make_model, sensor
import numpy as np
import time, math
def GetTb(Tz, H):
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
    return res.TbV()