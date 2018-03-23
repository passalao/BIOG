#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
from numpy import loadtxt

# Import SMOS data
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMean_TbV_52.5deg_xy.nc')
ny_Obs = Obs.dimensions['cols'].size
nx_Obs = Obs.dimensions['rows'].size
Tb_Obs = Obs.variables['BT_V']
Lon = Obs.variables['lon']
Lat = Obs.variables['lat']
Mask = Obs.variables['mask']

#Import GRISLI data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
H = GRISLI.variables['H']
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']

Sites=['DomeC', 'Vostok', 'DomeFuji', 'LawDome', 'T-UC', 'T-AIS']
Coords=np.array([[65,167],[86,165],[143,147],[62,206],[80,94],[79,93]])
Temp=np.zeros((np.shape(Sites)[0],5))

f=open("RT_Sites_"+str(BIOG.var.Perm)+'_'+str(BIOG.var.NbLayers)+"l.csv",'w')
f.write('Number of vertical layers: '+str(BIOG.var.NbLayers)+"\n")
f.write('Permittivity model: '+str(BIOG.var.Perm)+"\n")
f.write('Number of streams: '+str(BIOG.var.NbStreams)+"\n")
f.write('\n')
f.write('Site Tb Ts_LocalObs Ts_Dataset Tb_LocalObs Tb_Dataset \n')

i=0
for site in Sites:
    print(site)
    x=Coords[i,0]
    y=Coords[i,1]
    Tz_gr_at_Point = Tz_gr[x,y, :]
    Data = loadtxt("../../SourceData/Temperatures/"+str(site)+".csv", comments="#", delimiter=",",unpack=False)
    Tz_Obs = Data[:, 1]
    Z = Data[:, 0]
    Tb_modobs = BIOG.fun.GetTb(Tz_Obs[:], Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams,
                            BIOG.var.Perm, BIOG.var.RTModel)
    Tb_modgr = BIOG.fun.GetTb(Tz_gr_at_Point, Z[-1], BIOG.var.NbLayers, BIOG.var.Freq, BIOG.var.Angle, BIOG.var.NbStreams,
                           BIOG.var.Perm, BIOG.var.RTModel)
    Temp[i,0]= Tb_Obs[0,x,y]
    Temp[i,1]=Tz_Obs[0]+273.15
    Temp[i,2] =Tz_gr_at_Point[0]+273.15
    Temp[i,3] =Tb_modobs
    Temp[i,4] =Tb_modgr
    f.write(str(site)+' ')
    [f.write(str(Temp[i,j])+' ') for j in np.arange(0,5)]
    f.write("\n")
    i=i+1