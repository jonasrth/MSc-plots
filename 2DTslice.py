import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from scipy.interpolate import RegularGridInterpolator

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

iy = 100

lw = 1.0
alpha = 1.0

fs = np.array([12,15])*0.6
fig,ax = plt.subplots(5,4,figsize=(fs[0],fs[1]))

Teff_list = [6500,5770,5000,4000,3200]
B0_list = [0,100,200,500]

Tmax = []
Tmin = []
for i,Teff in enumerate(Teff_list):
    for j, B0 in enumerate(B0_list):
        F = h5py.File(input_files[Teff][B0], "r")
        T = np.array(F["temperature"][0,:,iy,:])
        Tmax.append(np.max(T))
        Tmin.append(np.min(T))
        F.close()

Tmax = np.log10(np.max(Tmax))
Tmin = np.log10(np.min(Tmin))

for i,Teff in enumerate(Teff_list):
    for j, B0 in enumerate(B0_list):
        F = h5py.File(input_files[Teff][B0], "r")
             
        T = np.array(F["temperature"][0,:,iy,:])
        #z = np.array(F["z"][0,0,0,:])*1e-8
        #x = np.cumsum(np.array(F["dx"]))*1e-8
        
         # Finding height of TR=1
        TR = np.array(F["tau"][0,:,iy,:])
        iTR1 = np.argmin(abs(TR-1), axis=1)
        iTRC = np.argmin(abs(TR-1e-5), axis=1)
        i_s = round(np.mean(iTR1))
        
        # Defining origin of z as height of TR=1
        z = np.array(F["z"][0,0,0,:])*1e-8
        z = z - z[i_s]
        x = np.cumsum(np.array(F["dx"]))*1e-8
        
        ax[i,j].plot(x, z[iTR1], color="white", ls="dashed", lw=lw, alpha=alpha)
        ax[i,j].plot(x, z[iTRC], color="white", ls="dashdot", lw=lw, alpha=alpha)
        
        if B0==0:
            x0 = 0.1*x[-1]
            y0 = 0.1*abs(z[0]-z[-1]) + z[-1]
            ax[i,j].plot([x0, x0+1], [y0,y0], lw=5, color="black")
            ax[i,j].plot([x0, x0+1], [y0,y0], lw=3, color="white")
        
        # Finding height of plasma beta unity
        if B0!=0:
            PG = np.array(F["Pgas"][0,:,iy,:])
            B2 = np.array(F["bx"][0,:,iy,:])**2 + np.array(F["by"][0,:,iy,:])**2 + np.array(F["bz"][0,:,iy,:])**2           
            PB = B2/(8*np.pi)
            beta = PG/PB
            i_beta = np.argmin(abs(beta-1), axis=1)
            #ax[i,j].plot(x, z[i_beta], color="white", ls="dotted", lw=lw, alpha=alpha)
        
        if Teff==3200:
            
            dz = abs(z[0]-z[1])
            DZ = abs(z[0]-z[-1])
            nz = int(DZ/dz)

            z_new = (np.arange(nz) - nz)*dz + z[0]
            
            X,Z = np.meshgrid(x,z_new, indexing="ij")
            
            interp = RegularGridInterpolator((x,z[::-1]), T[:,::-1])
            T = interp((X,Z), method="linear")[:,::-1]
            z = z_new[::-1]  
        
        logT = np.log10(T[:,::-1].T)
        
        ax[i,j].get_xaxis().set_ticks([])
        if j!=0:
            ax[i,j].get_yaxis().set_ticks([])
            
        extent = [x[0],x[-1], z[-1], z[0]]
        im = ax[i,j].imshow(logT, extent=extent, aspect="auto", origin="lower", vmin=Tmin, vmax=Tmax, cmap="inferno")
        F.close()
        
for i,Teff in enumerate(Teff_list):
    ax[i,-1].yaxis.set_label_position("right")
    ax[i,-1].set_ylabel(r"$T_{eff}$ = "+f"{Teff} K")
    ax[i,0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[i,0].set_ylabel(r"$z$ [Mm]")
for j, B0 in enumerate(B0_list):
    ax[-1,j].set_xlabel(fr"$B_0$ = {B0} G")

fig.colorbar(im, ax=ax.ravel().tolist(), location="top", fraction=0.046, pad=0.02, label=r"$\log_{10}(T$ [K]$)$")
plt.subplots_adjust(top=0.84, hspace=0, wspace=0)

figname = "2DTslice.pdf"
#plt.savefig("figures/"+figname, bbox_inches="tight")
plt.savefig("figures/presentation/"+figname, bbox_inches="tight")