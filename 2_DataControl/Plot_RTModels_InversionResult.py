###############################################################
#Plot the imaginary part of permittivities against temperature#
###############################################################
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math

##############################################################
#Models
def Permittivity(Model,T, D, Freq):
    e_ice=np.zeros((2,np.size(T)))

    if Model=="Tiuri":
        e_ice[0,:]=1+1.7*D/1000+0.7*(D/1000)**2
        e_ice[1,:]=1.5871e6*(0.52*D/1000+0.62*(D/1000)**2)*(1/Freq+1.23e-14*Freq**0.5)
        expT=np.exp(0.036*(T-273.15))
        e_ice[1,:]=e_ice[1,:]*expT

    if Model=="Matzler":
        e_ice[0,:] = 3.1884 + 9.1e-4 * (T - 273.0)
        theta = 300.0 / T - 1.0
        alpha = (0.00504 + 0.0062 * theta) * np.exp(-22.1 * theta)
        B1 = 0.0207
        B2 = 1.16e-11
        b = 335
        deltabeta = np.exp(-9.963 + 0.0372 * (T - 273.16))
        betam = (B1 / T) * (np.exp(b / T) / ((np.exp(b / T) - 1) ** 2)) + B2 * (Freq/1e9) ** 2
        beta = betam + deltabeta
        e_ice[1, :] = alpha /(Freq/1e9) + beta * (Freq/1e9)
    return e_ice

##############################################################
#Compute the permittivities, from direct model or inversion
c=299792458
freq=1.413e9
l=c/freq
rho=917

T=np.arange(210,270,1)
Im_ETiuri=Permittivity("Tiuri", T, rho, freq)[1]
Im_EMatzler=Permittivity("Matzler", T, rho, freq)[1]
Real_E=Permittivity("Tiuri", T, rho, freq)[0]

StartL=25/np.exp(0.036*(T-273.15))
StartE=l/(4*math.pi*StartL)*(Real_E[0])**0.5

#T_inv=np.array([-57.5,-50,-45,-40,-35,-30])+273.15+2.5#for East
#L_inv=np.array([448,301,331,286,255,258]) #East Antarctica

T_inv=np.arange(-60,-25,5)+2.5
L_inv=np.array([396,502,421,357,305,288,273]) #East Antarctica

E_inv=l/(4*math.pi*L_inv)*(Real_E[0])**0.5

T_inv2=np.arange(-45,-20,5)+2.5
print(T_inv2)
L_inv2=np.array([446,240,298,163,141]) #West Antarctica
E_inv2=l/(4*math.pi*L_inv2)*(Real_E[0])**0.5

#Discrepancy between retrieval and Mätzler
Im_EMatzler_controlEast=Permittivity("Matzler", T_inv, rho, freq)[1]
DiscE=(Im_EMatzler_controlEast-E_inv)#/np.size(T_inv)
Im_EMatzler_controlWest=Permittivity("Matzler", T_inv2, rho, freq)[1]
DiscW=(Im_EMatzler_controlWest[0:3]-E_inv2[0:3])#/np.size(E_inv2[0:2])
Disc=np.concatenate((DiscE, DiscW))
Disc=Disc/np.size(Disc)
Discrepancy=(np.dot(Disc, Disc))**0.5

print(Discrepancy)
##############################################################
#Plot
plt.plot(T-273.15, StartE*1e3, c='k', linestyle='--', linewidth=0.5, label="Initial guess")
plt.plot(T-273.15, Im_ETiuri*1e3, c='darkorange', label="Tiuri")
#plt.plot(T, Im_EMatzler*1e3*1.2, c='b', label="Mätzler")
plt.plot(T-273.15, Im_EMatzler*1e3, c='b', label="Mätzler")

plt.scatter(T_inv, E_inv*1e3, label="From retrieval, East", marker="D", s=20, c="r")
plt.scatter(T_inv2, E_inv2*1e3, label="From retrieval, West", marker="D", s=20, c="darkgreen")
plt.xlabel("Ice temperature (K)", fontsize=15)
plt.ylabel(r"$\epsilon'' \times$ 1e3", fontsize=15)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlim([-60,-20])
plt.ylim([0,1.2])
plt.legend(fontsize=13)
plt.grid()
plt.show()