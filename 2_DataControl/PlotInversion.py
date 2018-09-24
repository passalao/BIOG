#!/usr/bin/python
# -*- coding: cp1252 -*-

###############################################################
#Plot the evolution of the cost function                      #
###############################################################
import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt

##############################################################
#Import data  of retrieval process
Data = loadtxt("../../OutputData/InversionProcess_40-35East.csv", comments="#", delimiter=" ", unpack=False)
Depth=Data[:,0]
J=(Data[:,1])**0.5
R=(Data[:,2])**0.5
E=Data[:,3]
iterate=np.arange(0, np.size(Depth), 1)

##############################################################
#Plot

f, (ax1, ax2, ax3, ax4) = plt.subplots(1,4, sharex=True)

ax1.plot(iterate, J, 'coral', linewidth=0.5, label='$\sqrt{\mathcal{J}}$ (K)' )
ax1.set_ylim([0,2])
ax1.legend(loc=3)
ax1.grid()
ax1.set_title("a)")

ax2.plot(iterate, R, 'g', linewidth=0.5, label='$\sqrt{\mathcal{R}}$')
ax2.set_ylim([0,0.5])
ax2.set_xlim([0,np.size(iterate)-1])
ax2.legend(loc=3)
ax2.grid()
ax2.set_title("b)")
ax2.set_xlabel('Iteration')

ax3.plot(iterate, E, 'purple', linewidth=0.5, label='$\overline{\eta}$')
ax3.set_ylim([0.97,0.98])
ax3.legend(loc=3)
ax3.grid()
ax3.set_title("c)")

#ax4 and ax5 are both on the same plot
ax4.plot(iterate, 1/Depth, 'b', linewidth=0.5, label='$\kappa_a\, (m^{-1})$')
ax4.grid()
ax4.set_title("d)")
ax5=ax4.twinx()
ax5.plot(iterate, Depth, 'r', linewidth=0.5, label='$1/\kappa_a\, (m)$')
ax4.legend(loc=3)
ax5.legend(loc=2)

#Rotate ticks to save space on the plot
plt.setp(ax1.get_yticklabels(), rotation = 90)
plt.setp(ax2.get_yticklabels(), rotation = 90)
plt.setp(ax3.get_yticklabels(), rotation = 90)
plt.setp(ax4.get_yticklabels(), rotation = 90)
plt.setp(ax5.get_yticklabels(), rotation = 90)

plt.show()