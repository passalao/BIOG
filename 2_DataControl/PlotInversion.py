#!/usr/bin/python
# -*- coding: cp1252 -*-

###############################################################
#Plot the cost function evoluion                              #
###############################################################
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import loadtxt

##############################################################
#Import data  of inversion process
Data = loadtxt("../../OutputData/InversionProcess_40-35East.csv", comments="#", delimiter=" ", unpack=False)
Depth=Data[:,0]
J1=(Data[:,1])**0.5
J3=(Data[:,2])**0.5
E=Data[:,3]
iterate=np.arange(0, np.size(Depth), 1)

##############################################################
#Plot
#plt.plot(iterate, Depth, linewidth=1)

f, (ax2, ax3, ax5, ax1) = plt.subplots(1,4, sharex=True)
ax1.plot(iterate, 1/Depth, 'b', linewidth=0.5, label='$\kappa_a\, (m^{-1})$')
ax1.grid()
ax1.set_title("d)")

ax4=ax1.twinx()
ax4.plot(iterate, Depth, 'r', linewidth=0.5, label='$1/\kappa_a\, (m)$')
ax1.legend(loc=3)
ax4.legend(loc=2)


ax5.plot(iterate, E, 'purple', linewidth=0.5, label='$\overline{\eta}$')
#ax5.set_ylim([0.95,1])
ax5.set_ylim([0.97,0.98])
ax5.legend(loc=3)
ax5.grid()
ax5.set_title("c)")

ax2.plot(iterate, J1, 'coral', linewidth=0.5, label='$\sqrt{\mathcal{J}}$ (K)' )
ax3.set_xlabel('Iteration')
ax2.set_ylim([0,2])
ax2.legend(loc=3)
ax2.grid()
ax2.set_title("a)")

ax3.plot(iterate, J3, 'g', linewidth=0.5, label='$\sqrt{\mathcal{R}}$')
ax3.set_ylim([0,0.5])
ax3.set_xlim([0,np.size(iterate)-1])
ax3.legend(loc=3)
ax3.grid()
ax3.set_title("b)")



plt.setp(ax1.get_yticklabels(), rotation = 90)
plt.setp(ax2.get_yticklabels(), rotation = 90)
plt.setp(ax3.get_yticklabels(), rotation = 90)
plt.setp(ax4.get_yticklabels(), rotation = 90)
plt.setp(ax5.get_yticklabels(), rotation = 90)

#plt.yticks(rotation='vertical')

plt.show()