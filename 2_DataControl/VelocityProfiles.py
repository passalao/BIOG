#!/usr/bin/python
# -*- coding: cp1252 -*-
#
#Ce script lit un fichier csv
#et trace les valeurs
#

#Bibliotheques
import os, csv, math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
axis_font = { 'size':'18'}

zeta = np.linspace(0.,1., 101)
p3=5.
omega3=np.zeros(len(zeta))

i=0
for z in zeta:
	omega3[i]=1-(p3+2)/(p3+1)*(1-zeta[i])+1/(p3+1)*(1-zeta[i])**(p3+2)
	i=i+1
H1=3500
H2=2500
Z1=(1-zeta)*H1
Z2=(1-zeta)*H2
#######################################################################
#Affichage du graphique

fig = plt.figure()
plt.grid()

#plt.plot(omega1,1-zeta,'k', linewidth=0.9, label="$p=1$")
#plt.plot(omega2,1-zeta,'k--', linewidth=0.9, label="$p=2.5$")
plt.plot(omega3,Z1,'b-.', linewidth=0.9, label="$p=5$")
plt.plot(omega3,Z2,'r-.', linewidth=0.9, label="$p=5$")
plt.plot([1,0],[0,H1],c='b')
plt.plot([1,0],[0,H2],c='r')
#plt.plot(omega4,1-zeta,'k:', linewidth=0.9, label="$p=10$")
plt.legend(loc=2)
plt.tick_params(labelsize=14)
plt.xlabel('$\omega$ - $\mathrm{Shape\,function}$', fontsize=16)
plt.ylabel("$\zeta $- $\mathrm{Normalized\,depth}$", fontsize=16)
#plt.xlim(0.0,1.0)
#plt.ylim(0.0,1.0)
ax = plt.gca()
ax.invert_yaxis()
plt.show()

