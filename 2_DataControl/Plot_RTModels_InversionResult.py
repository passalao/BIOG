import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math

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

freq=1.413e9
l=299792458/freq
rho=917

T=np.arange(210,270,1)
Im_ETiuri=Permittivity("Tiuri", T, rho, freq)[1]
Im_EMatzler=Permittivity("Matzler", T, rho, freq)[1]
Real_E=Permittivity("Tiuri", T, rho, freq)[0]

#T_inv=np.arange(-60,-10,5)+273.15+2.5
#L_inv=np.array([633,412,332,282,241,206,177,152,132,115])
T_inv=np.arange(-60,-15,5)+273.15+2.5
L_inv=np.array([800,770,486,423,435,356,236,192,163])#633,412,332,282,241,206,177,152,132,115])

E_inv=1/L_inv
E_inv=l/(4*math.pi*L_inv)*(Real_E[0])**0.5

#L=l/(4*math.pi)/Im_EMatzler*(Real_E[0])**0.5

#E_inv=5e-8*np.exp(0.036*T_inv)

plt.plot(T, Im_ETiuri*1e3, c='r', label="Tiuri")
plt.plot(T, Im_EMatzler*1e3, c='b', label="Mätzler")
plt.scatter(T_inv, E_inv*1e3, label="From inversion", marker="+", s=50, c="darkorange")
plt.xlabel("Ice temperature (K)")
plt.ylabel(r"$\epsilon'' \times$ 1e3")
plt.legend()
plt.grid()
plt.show()