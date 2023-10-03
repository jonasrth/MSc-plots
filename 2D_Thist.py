import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

#HP_dict = np.load("files/HP_dict.pkl", allow_pickle=True)
TR_dict = np.load("files/TR_dict.pkl", allow_pickle=True)

def plot_T_hist(Teff, Tbins=200):
    """
    """
    
    fig,ax = plt.subplots(2, figsize=(10,8), sharex=True)
    
    letter_list = ["a", "b", "c"]
    
    zmax = 100
    for i,B0 in enumerate([0, 500]):
        
        F = h5py.File(input_files[Teff][B0], "r")
    
        T = np.array(F["temperature"][0,...])
        z = np.array(F["z"][0,...])*1e-8
        
        logTR = TR_dict[Teff][B0]
        
        i_s = np.argmin(abs(logTR))
        i_c = np.argmin(abs(logTR - (-5)))
        
        if B0==0:
            z_ = z - z[0,0,i_s]
            
        ax[i].axvline(x = z_[0,0,i_s], color="k", ls="dashed", lw=1.2, zorder=3)
        ax[i].axvline(x = z_[0,0,i_c], color="k", ls="dashdot", lw=1.2, zorder=3)
        
        if z[0,0,0] < zmax:
            zmax = z[0,0,0]
    
        if Teff==3200:
            """
            Grid cells for Teff=3200 K model are not uniform.
            """

            dx = np.array(F["dx"])*1e-8
            x = np.cumsum(dx)
            nx = len(x)

            dz = abs(z[0,0,0]-z[0,0,1])
            DZ = abs(z[0,0,0]-z[0,0,-1])
            nz = int(DZ/dz)

            z_new = (np.arange(nz) - nz)*dz + z[0,0,0]

            X,Y,Z = np.meshgrid(x,x,z_new)

            interp = RegularGridInterpolator((x,x,z[0,0,::-1]), T[...,::-1])
            T = interp((X,Y,Z))
            z = Z
    
        zbins = z.shape[-1]
        H, x_edges, y_edges = np.histogram2d(z.flatten(), T.flatten(), bins=[zbins,Tbins])
    
        ax[i].grid(alpha=1, zorder=1)

        extent = [z_[0,0,-1], z_[0,0,0], np.min(T)*1e-3, np.max(T)*1e-3]
        #extent = [0, 1, np.min(T)*1e-3, np.max(T)*1e-3]
        ax[i].imshow(np.log10(H.T+1), extent=extent, origin="lower", aspect="auto", cmap="Greys", alpha=0.75, zorder=2)
        ax[i].plot(z_[0,0,:], np.mean(T*1e-3, axis=(0,1)), color="red", ls="dashed", zorder=3)
        ax[i].set_ylabel("Gas temperature [kK]")
        ax[i].text(0.9, 0.8, letter_list[i], color="black", size=20, transform=ax[i].transAxes)
        
    for i in range(2):
        #ax[i].set_xlim([z_[0,0,-1], zmax])
        ax[i].set_xlim([z_[0,0,-1], z_[0,0,0]])
    ax[-1].set_xlabel("z [Mm]")
    
    plt.subplots_adjust(hspace=0.1)
    
    figname = "2D_Thist.pdf"
    plt.savefig("figures/"+figname, bbox_inches="tight")
    
plot_T_hist(5000)
