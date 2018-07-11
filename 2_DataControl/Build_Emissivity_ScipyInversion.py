from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, time
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG
import NC_Resources as ncr
import scipy.optimize as opt

def InitEmissivity(Depth, frac):#Tz_gr,H,Depths,f,TsMin, TsMax,Subsample):
    Teff=np.zeros(np.shape(H))#+1e10
    TbObs=np.zeros(np.shape(H))
    Mask=np.zeros(np.shape(H))

    for i in np.arange(0,np.shape(H)[0], 1):
        for j in np.arange(0,np.shape(H)[1], 1):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>10:
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)
                Mask[i,j]=1
    TbObs=Tb*Mask
    Emissivity=TbObs/Teff
    #Emissivity[Teff==0]=0
    #print(Emissivity)
    Rand=np.random.normal(0, 0.001, np.size(Tb))
    Rand=np.reshape(Rand,(np.shape(Tb)))
    Ebarre=np.mean(Emissivity[Emissivity!=0])
    Emissivity=frac*Ebarre+(1-frac)*Emissivity+Rand
    return Emissivity

def ComputeLagComponents(x):
    #print(x)
    Depth=x[0]
    Mu=x[1]
    E=x[2:]
    E=np.reshape(E,np.shape(Tb))

    Teff=np.zeros(np.shape(Tb))
    Mask=np.zeros(np.shape(Tb))

    for i in np.arange(0,np.shape(H)[0], 1):
        for j in np.arange(0,np.shape(H)[1], 1):
            if Ts[i,j]>TsMin and Ts[i,j]<TsMax and H[i,j]>1:
                Teff[i,j]=273.15+sum(Tz_gr[i,j]*np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)/sum(np.exp(-(1-Zeta[i,j])*H[i,j]/Depth)*0.05*H[i,j]/Depth)
                Mask[i,j]=1
    TbObs=Tb*Mask
    #Tbmod=np.reshape(E,np.shape(Tb))*Teff
    Tbmod=E*Teff

    Teff1D=np.reshape(Teff, (1,np.size(Teff)))[0,:]
    Tbmod1D=np.reshape(Tbmod, (1,np.size(Tbmod)))[0,:]
    TbObs1D=np.reshape(TbObs, (1,np.size(TbObs)))[0,:]
    Emiss1D=np.reshape(E, (1,np.size(Emissivity)))[0,:]

    J1=np.sum((Tbmod1D[Tbmod1D!=0]-TbObs1D[Tbmod1D!=0])**2)/np.size(Tbmod1D)#np.std(TbObs1D[Tbmod1D!=0])**2
    #J2=np.std(Emiss1D[Emiss1D!=0])

    #Compute normalized covariance = correlation
    m=np.stack((Emiss1D[Emiss1D!=0], Teff1D[Emiss1D!=0]), axis=0)
    J3=abs(np.cov(m)[0,1])/np.std(Emiss1D[Emiss1D!=0])/np.std(Teff1D[Emiss1D!=0])#)**2)**0.5
    print(J1+Mu*J3)
    return J1+Mu*J3#, J2, J3, Teff, Tbmod, TbObs, Teff1D, Tbmod1D, TbObs1D, Emiss1D, Mask

# Import SMOS data
print("Load data")
Obs = netCDF4.Dataset('../../SourceData/SMOS/SMOSL3_StereoPolar_AnnualMeanSansNDJ_TbV_52.5deg_xy.nc')
Tb = Obs.variables['BT_V']
X = Obs.variables['x_ease2']
Y = Obs.variables['y_ease2']
Tb=Tb[0]
nc_obsattrs, nc_obsdims, nc_obsvars = ncr.ncdump(Obs)

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4Ts.nc')
H = np.array(GRISLI.variables['H'])
Zeta = GRISLI.variables['Zeta']
Tz_gr = GRISLI.variables['T']
Ts=Tz_gr[:,:,0]

#Parameters for data analyse
DeltaT=5
SliceT=np.arange(-60, -20,DeltaT)
print("Temperatures : ", SliceT)

FinalEmissivity=np.zeros(np.shape(Ts))
Depths=[]
J1s=[]
J3s=[]

#Work slice by slice
for t in SliceT:
    TsMax = t + DeltaT
    TsMin = t
    print("  ")
    print("Slice between ", TsMin, " and ", TsMax)

    Depth=-5*t #initiate with plausible depth
    Mu=50
    Emissivity=np.random.normal(0.98, 0.02, np.size(Ts))
    Emissivity=np.reshape(Emissivity,(np.shape(Ts)))
    #Emissivity=InitEmissivity(Depth, 0.5)

    Emissivity=np.reshape(Emissivity,(1,np.size(Emissivity)))
    print(np.shape([[Depth]]), np.shape(Emissivity[0]))
    x0=np.concatenate(([Depth],[Mu], Emissivity[0]), axis=0)
    epsilon=np.concatenate(([10], [100], np.ones(np.size(Emissivity))*5e-3), axis=0)
    print(x0)

    #here optimization with scipy
    BestValue=opt.fmin(ComputeLagComponents,x0)#, method='Newton-CG')#, fprime=None)#, epsfcn=epsilon)
    print(BestValue)

    print("Lambda", Lambda)
    print("Mu:", Mu)
    Depths.append(Depth)
    J1s.append(J1)
    J3s.append(J3)

    FinalEmissivity=FinalEmissivity+Mask*Emissivity

print("Depths:" ,Depths)
print("J1:", J1s)
print("J3:", J3s)

#for display
FinalEmissivity[FinalEmissivity==0]=FinalEmissivity[112,100]

#Create NetCDF file
outfile = r'../../SourceData/WorkingFiles/Emissivity_FromGradientDescent_nDim.nc'
nc_new = netCDF4.Dataset(outfile, 'w', clobber=True)

cols = len(X[0,:])
rows = len(Y[:,0])

Xout = nc_new.createDimension('x', cols)
Xout = nc_new.createVariable('x', 'f4', ('x',))
Xout.standard_name = 'x'
Xout.units = 'm'
Xout.axis = "X"
Xout[:] = X[0,:]-25000
Yout = nc_new.createDimension('y', rows)
Yout = nc_new.createVariable('y', 'f4', ('y',))
Yout.standard_name = 'y'
Yout.units = 'm'
Yout.axis = "Y"
Yout[:] = Y[:,0]+25000

nc_new.createVariable("Emissivity", 'float64', ('y','x'))
nc_new.variables["Emissivity"][:] = FinalEmissivity[::-1,:]
crs = nc_new.createVariable('spatial_ref', 'i4')
crs.spatial_ref='PROJCS["WGS_84_NSIDC_EASE_Grid_2_0_South",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",-90],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
nc_new.close()

fig, ax = plt.subplots(nrows=1, ncols=1)
norm = mpl.colors.Normalize(vmin=0.9, vmax=1)
cmap = mpl.cm.spectral
myplot = ax.pcolormesh(FinalEmissivity, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=np.arange(0.90, 1.01, 0.02))
cbar.set_label('Emissivity', rotation=270)
cbar.ax.set_xticklabels(['0.95', '0.96', '0.97', '0.98', '0.99', '1.0'])
plt.savefig("../../OutputData/img/InvertingEmissDepth/Emissivity_DescentGrad_nDim.png")
plt.show()