from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyproj
from sklearn import linear_model
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

def ComputeCorr(Tz_gr,H,Depth,TsMin, TsMax,Subsample):
    Teff = []
    TbObs=[]
    for i in np.arange(0,np.shape(H)[0], Subsample):
        for j in np.arange(0,np.shape(H)[1], Subsample):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff.append(273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j,:])*H[i,j]/Depth)*0.05*H[i,j]/Depth))
                TbObs.append(Tb[i,j])
    Teff=np.array(Teff)
    TbObs=np.array(TbObs)
    Teff=Teff[TbObs < 1e45]
    TbObs = TbObs[TbObs < 1e45]
    Emissivity=np.array(TbObs)/np.array(Teff)
    emoy = np.mean(Emissivity)
    J1=sum((Teff-TbObs)**2)
    J2=sum((Emissivity[Emissivity>=1.0]-emoy)**2)
    #print("Depth:", Depth, "Emiss:", np.mean(Emissivity))

    #Compyute regression
    regr = linear_model.LinearRegression()
    regr.fit(TbObs[:, np.newaxis], Teff)
    x_test = np.linspace(np.min(TbObs), np.max(TbObs), 100)
    print("Score:", regr.score(TbObs[:, np.newaxis], Teff))

    #if Depth==700:
    '''plt.plot(x_test, regr.predict(x_test[:, np.newaxis]), color='blue')
    plt.scatter(TbObs, Teff,s=0.1)
    plt.plot([210,270],[210,270],c='r')
    plt.text(220,260,d)
    plt.show()'''
    return Teff, J1, J2

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Mask = Obs.variables['mask']
Tb=Tb[0]

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4Ts.nc')
H = np.array(GRISLI.variables['H'])
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

TsMax=-45
TsMin=-50
Subsample=1
LayerThick=10

#Determine which Depth is the best one
Depths=np.arange(300,1100,100)
#Alpha=np.arange(0,-1e6,-1e4)

#for a in Alpha:
J1 = []
J2 = []
for d in Depths:
    Data=ComputeCorr(Tz_gr, H, d, TsMin, TsMax, Subsample)
    J1.append(Data[1])
    J2.append(Data[2])
    #print((max(J1)-min(J1))/(max(J2)-min(J2)))
    #Jtot=J1[-1]+a*J2[-1]
    #print(a,d,Jtot)
print(min(J1), max(J1))
J1=np.array(J1)
J2=np.array(J2)
J1=(J1-min(J1))/(max(J1)-min(J1))
J2=(J2-min(J2))/(max(J2)-min(J2))

Jtot=(J1+J2)/2

plt.plot(Depths,J1,c="b", label="J1: Teff-Tb")
plt.plot(Depths,J2,c="r", label="J2: Emiss-Emissbar")
plt.plot(Depths,Jtot,c="k", label='Total')
plt.xlabel('Depths (m)')
plt.ylabel('Cost functions, normalized')
plt.legend()
plt.title(str(TsMin)+ "< Ts < "+str(TsMax))
plt.show()

print("Compute Tmoy")


TmoyTb=Tmoy-Tb

#Calculate dependence of Ts-Tb on Ts-TatDepth
Tb1D=np.reshape(Tb,(np.size(DeltaT),1))
Tmoy1D=np.reshape(Tmoy,(np.size(DeltaT),1))
DeltaT1D=np.reshape(DeltaT,(np.size(DeltaT),1))
TbTs1D=np.reshape(TbTs,(np.size(Tb),1))
TmoyTb1D=np.reshape(TmoyTb,(np.size(Tb),1))

Ts1D=np.reshape(TsRef,(np.size(Tb),1))
Tb1D=np.reshape(Tb,(np.size(Tb),1))
Mask1D=np.reshape(Mask,(np.size(Tb),1))

#Clean dataset from useless points
i=0
DeltaT1D_=[]
TbTs1D_=[]
TmoyTb1D_=[]
for d in DeltaT1D:
   tm=TmoyTb1D[i]
   tbts=TbTs1D[i]
   ts=Ts1D[i]
   if d >-50 and tm<220 and ts>-35 and ts <-30:
       DeltaT1D_.append(d)
       TmoyTb1D_.append(tm)
       TbTs1D_.append(tbts)
   i=i+1

#Least square resolution
# We're solving Ax = B
A = np.column_stack([np.ones(len(DeltaT1D_)), DeltaT1D_])
B = TbTs1D_

# Solve the system of equations.
result,_,_,_= np.linalg.lstsq(A,B, rcond=None)
a, b = result
print(result)

#Correction for a
a=0

#Compute the glaciologic and radiative component of Tb-Ts
TbTsTrend1D=b*(DeltaT1D)
TbTsModul1D=TbTs1D-TbTsTrend1D
#Reshape data, and correct for the difference between the Ts datasets
TbTsModul=np.reshape(TbTsModul1D,np.shape(TatDepth))+b*(TsRef-TsRACMO)
TbTsTrend=np.reshape(TbTsTrend1D,np.shape(TatDepth))-b*(TsRef-TsRACMO)
Emissivity=(TbTsModul+TbTsTrend+(TsRef+273.15))/(Tmoy+TsRef-TsRACMO)#Calculé dans l'univers "Crocus"
Emissivity2=Tb/(Tmoy+TsRef-TsRACMO)#Identical to Emissivity. So the results are consistent whatever the way emissivity is computed.

Em1D=np.reshape(Emissivity, (np.size(Tb),1))

TbTsModul=TbTsModul*(4-np.array(Mask))/3
TbTsTrend=TbTsTrend*(4-np.array(Mask))/3
Emissivity=Emissivity*(4-np.array(Mask))/3

#Scatter
cmap = plt.cm.get_cmap("coolwarm")
norm=mpl.colors.Normalize(vmin=0.95, vmax=1.05)
plt.scatter(DeltaT1D,TbTs1D, s=5, c=Em1D, cmap=cmap, norm=norm, edgecolors='none')
plt.plot([-10,20],[a+b*(-10),a+20*b],color='r', lw=0.5)
plt.xlabel("Tmoy - Ts (K)")
plt.ylabel("Tb - Ts (K)")
plt.xlim(-10,10)
plt.ylim(-20,20)
#plt.colorbar()
cbar = plt.colorbar()
cbar.set_label('Emissivity', rotation=270)
plt.grid()
#plt.show()

'''# Geographic plot
fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)
cmap = mpl.cm.spectral
norm = mpl.colors.Normalize(vmin=-15, vmax=15)
myplot=ax[0].pcolormesh(TbTs, cmap=cmap, norm=norm)
ax[0].set_title('Tb-Ts')
ax[1].pcolormesh(TbTsTrend, cmap=cmap, norm=norm)
ax[1].set_title('Regression : Tmoy-Ts')
ax[2].pcolormesh(TbTsModul, cmap=cmap, norm=norm)
ax[2].set_title('Obs-Reg : Tb-Tmoy')
plt.axis('equal')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
cbar = fig.colorbar(myplot, cax=cbar_ax, ticks=np.arange(-15, 16, 5))
cbar.set_label('Kelvin difference (K)', rotation=270)
cbar.ax.set_xticklabels(['-15', '-5', '5', '15'])
plt.show()'''

'''fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0.95, vmax=1)
cmap = mpl.cm.spectral
myplot=ax.pcolormesh(Emissivity, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0.95, 1.01, 0.01))
cbar.set_label('Emissivity', rotation=270)
cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
plt.show()'''

'''#Write output NetCDF file
outfile = r'../../SourceData/WorkingFiles/TbDeconvolution.nc'
cols = np.shape(H)[1]
rows = np.shape(H)[0]
dsout = netCDF4.Dataset(outfile, 'w', clobber=True)
Yout = dsout.createDimension('y', rows)
Yout = dsout.createVariable('y', 'f4', ('y',))
Yout.standard_name = 'y'
Yout.units = 'm'
Yout.axis = "Y"
Yout[:] = Y[:,0]+25000
Xout = dsout.createDimension('x', cols)
Xout = dsout.createVariable('x', 'f4', ('x',))
Xout.standard_name = 'x'
Xout.units = 'm'
Xout.axis = "X"
Xout[:] = X[0,:]-25000

dsout.createVariable('Ref-RACMO','float64',('y','x'))
dsout.variables['Ref-RACMO'][:] = np.array(TsRef[::-1,:]-TsRACMO[::-1,:])
dsout.createVariable('Ts-Tb','float64',('y','x'))
dsout.variables['Ts-Tb'][:] = np.array(TbTs[::-1,:])
dsout.createVariable('TbTsTrend','float64',('y','x'))
dsout.variables['TbTsTrend'][:] = np.array(TbTsTrend[::-1,:])
dsout.createVariable('TbTsModul','float64',('y','x'))
dsout.variables['TbTsModul'][:] = np.array(TbTsModul[::-1,:])
dsout.createVariable('Emissivity','float64',('y','x'))
dsout.variables['Emissivity'][:] = np.array(Emissivity[::-1,:])
crs = dsout.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
'''