#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

#Fictitious point, for verification
print("         Fictitious point         ")
T=[-53.15,-52.15]
H=3000.

angles=np.linspace(40,65,21)
Tb_gr=np.zeros((21,2))

i=0
for t in T:
    #Tz_gr_at_Point=t*np.ones(10)
    Tz_gr_at_Point = np.linspace(t,t,10)
    for a in angles:
        Tb_gr[np.where(angles==a)[0],i] = BIOG.fun.GetTb(Tz_gr_at_Point, H, BIOG.var.NbLayers, BIOG.var.Freq, a, BIOG.var.NbStreams, "Tiuri", "DMRT-ML")
    plt.plot(angles, Tb_gr[:,i])
    i=i+1
#print("Brewster angle:", np.where(Tb_gr==max(Tb_gr))[0])
print(Tb_gr[:,:])
plt.show()
plt.close()
