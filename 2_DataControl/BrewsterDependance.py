#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import matplotlib.pyplot as plt
import numpy as np
import sys, math
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

#Fictitious point, for verification
print("         Fictitious point         ")
T=[-53.15,-52.15]
H=3000.

angles=np.linspace(40,65,21)
Tb_gr=np.zeros((21,2))
NbLayers=[10]#,30, 100, 300]
NbStreams=[8,16,32,64]#,128]

'''for n in NbLayers:
    Tz_gr_at_Point = np.zeros(n)
    print("n=",n)
    for a in angles:
        depth=np.linspace(0,H,n)
        for d in depth:
            i=np.where(depth==d)[0]
            Tz_gr_at_Point[i] = -2-T[0] * math.exp(-(H - d) / (H / 6)) + T[0]
        #plt.plot(Tz_gr_at_Point,depth)
        #plt.show()
        #Tz_gr_at_Point = np.linspace(T[0], T[0]+20, n)
        Tb_gr[np.where(angles==a)[0],0] = BIOG.fun.GetTb(Tz_gr_at_Point, H, n, BIOG.var.Freq, a, BIOG.var.NbStreams, "Tiuri", "DMRT-ML")
        for d in depth:
            i=np.where(depth==d)[0]
            Tz_gr_at_Point[i] = -2-T[1] * math.exp(-(H - d) / (H / 6)) + T[1]
        #Tz_gr_at_Point = np.linspace(T[1], T[1]+20, n)
        Tb_gr[np.where(angles==a)[0],1] = BIOG.fun.GetTb(Tz_gr_at_Point, H, n, BIOG.var.Freq, a, BIOG.var.NbStreams, "Tiuri", "DMRT-ML")
    plt.plot(angles, Tb_gr[:,1]-Tb_gr[:,0], label=str(n)+' layers')'''

for s in NbStreams:
    for n in NbLayers:
        Tz_gr_at_Point = np.zeros(n)
        print("n=", n)
        for a in angles:
            depth = np.linspace(0, H, n)
            for d in depth:
                i = np.where(depth == d)[0]
                Tz_gr_at_Point[i] = -2 - T[0] * math.exp(-(H - d) / (H / 6)) + T[0]
            # plt.plot(Tz_gr_at_Point,depth)
            # plt.show()
            # Tz_gr_at_Point = np.linspace(T[0], T[0]+20, n)
            Tb_gr[np.where(angles == a)[0], 0] = BIOG.fun.GetTb(Tz_gr_at_Point, H, n, BIOG.var.Freq, a,
                                                                s, "Tiuri", "DMRT-ML")
            for d in depth:
                i = np.where(depth == d)[0]
                Tz_gr_at_Point[i] = -2 - T[1] * math.exp(-(H - d) / (H / 6)) + T[1]
            # Tz_gr_at_Point = np.linspace(T[1], T[1]+20, n)
            Tb_gr[np.where(angles == a)[0], 1] = BIOG.fun.GetTb(Tz_gr_at_Point, H, n, BIOG.var.Freq, a,
                                                                s, "Tiuri", "DMRT-ML")
        plt.plot(angles, Tb_gr[:,1]-Tb_gr[:,0], label=str(s)+' streams')

#plt.plot(angles, Tb_gr[:,:])
#plt.plot(angles, Tb_gr[:,1]-Tb_gr[:,0])
#print("Brewster angle:", np.where(Tb_gr==max(Tb_gr))[0])
plt.grid('on')
plt.xlabel('Incidence angle (deg)')
plt.ylabel('Emissivity')
plt.legend()
plt.show()
plt.close()
