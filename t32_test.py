import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator

from matplotlib.ticker import FormatStrFormatter

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

letter_list = ["a", "b", "c", "d"]
#TR_dict = np.load("../../Python/files/TR_dict.pkl", allow_pickle=True)

Teff = 3200
ix = 60

F = h5py.File(input_files[Teff][0], "r")

TR = np.array(F["tau"][0,ix,:,:])
iTR1 = np.argmin(abs(TR-1), axis=1)
i_s = round(np.mean(iTR1))

fs = np.array([10,10])*0.8
fig,ax = plt.subplots(2,2,figsize=(fs[0],fs[1]), sharey=True)

for i,B0 in enumerate([0,100,200,500]):
    F = h5py.File(input_files[Teff][B0], "r")
    
    z = np.array(F["z"][0,0,0,:])*1e-8
    x = np.cumsum(np.array(F["dx"]))*1e-8
    
    V = np.sqrt(np.array(F["vx"][0,ix,:,:])**2 + np.array(F["vy"][0,ix,:,:])**2 + np.array(F["vz"][0,ix,:,:])**2)*1e-5
    
    #TR = np.array(F["tau"][0,ix,:,:])
    #iTR1 = np.argmin(abs(TR-1), axis=1)
    iTRC = np.argmin(abs(TR-1e-5), axis=1)
    #i_s = round(np.mean(iTR1))
    
    z = z - z[i_s]
    
    
    if Teff==3200:
            
            dz = abs(z[0]-z[1])
            DZ = abs(z[0]-z[-1])
            nz = int(DZ/dz)

            z_new = (np.arange(nz) - nz)*dz + z[0]
            
            X,Z = np.meshgrid(x,z_new, indexing="ij")
            
            interp = RegularGridInterpolator((x,z[::-1]), V[:,::-1])
            V = interp((X,Z), method="linear")[:,::-1]
            z = z_new[::-1]
            
    logV = np.log10(V[:,::-1].T)
    
    extent = [x[0],x[-1],z[i_s],z[0]]
    im1 = ax.flat[i].imshow(logV[i_s:,:], origin="lower", aspect="auto", extent=extent,cmap="Reds")
    fig.colorbar(im1, ax=ax.flat[i], location="top", label=r"$\log_{10}$(v [km/s])")
    ax.flat[i].set_xlabel("y [Mm]")
    
    ax.flat[i].plot(x, z[iTRC], ls="dashdot", color="white")
    ax.flat[i].text(0.8,0.8,letter_list[i], color="white", size=20, transform=ax.flat[i].transAxes)
    #fig.colorbar(im2, ax=ax[i], location="top")

ax[0,0].set_ylabel("z [Mm]")
ax[1,0].set_ylabel("z [Mm]")
plt.subplots_adjust(wspace=0.05, hspace=0.25)

figname="t32_test.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")