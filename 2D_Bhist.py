import numpy as np
import matplotlib.pyplot as plt

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

def plot_2D_Bhist(Teff, B):
    """
    """
    
    fig,ax = plt.subplots(3, figsize=(10,12), sharex=True)
    
    F = h5py.File(input_files[Teff][B], "r")
    f = h5py.File(output_files[Teff][B][0.0005], "r")
    
    ### Adjust z origin to Rosseland mean opacity unity
    z = np.array(F["z"][0])*1e-8
    z_top = np.max(np.array(f["contfctz"]))*1e-8
    z = z - z[0,0,0] + z_top
    
    b_dict = {}
    b_dict["bz"] = abs(np.array(F["bz"][0]))
    b_dict["bxy"] = np.sqrt(np.array(F["bx"][0])**2 + np.array(F["by"][0])**2)
    b_dict["b"] = np.sqrt(np.array(F["bx"][0])**2 + np.array(F["by"][0])**2 + np.array(F["bz"][0])**2)
    
    ### Height of plasma beta unity
    Pb = np.mean(b_dict["b"]**2/(8*np.pi), axis=(0,1))
    Pg = np.mean(np.array(F["Pgas"][0]), axis=(0,1))
    beta = Pg/Pb
    zb = z[0,0,np.argmin(abs(beta - 1))]
    
    ### Histogram:
    
    ylabel_dict = {"bz"  : "Vertical $B$ [kG]", 
                   "bxy" : "Horisontal $B$ [kG]", 
                   "b"   : "Total $B$ [kG]"}
    letter_list = ["a", "b", "c"]
    for i,b in enumerate(b_dict.keys()):
        
        y = b_dict[b]
    
        yrange = (round(np.min(y/10))*10, round(np.max(y/10))*10)
        ybin_width = 10
        ybins = int((yrange[-1] - yrange[0])/ybin_width)
    
        H, xedges, yedges = np.histogram2d(z.flatten(), 
                                           y.flatten(),
                                           range=[[z[0,0,-1], z[0,0,0]], yrange],
                                           bins=[z.shape[-1], ybins])
        
        ax[i].grid(zorder=1)
        ax[i].imshow(np.log10(H.T+1), origin="lower", extent=[z[0,0,-1], z[0,0,0], yrange[0]*1e-3, yrange[-1]*1e-3], 
                                      aspect="auto", cmap="Greys", alpha=0.75, zorder=2)
        ax[i].plot(z[0,0,:], np.mean(y*1e-3, axis=(0,1)), color="red", ls="dashed", zorder=3)
        ax[i].axvline(x = zb, color="red", ls="dotted", zorder=3)
        ax[i].set_ylabel(ylabel_dict[b])
        ax[i].text(0.9, 0.8, letter_list[i], color="black", size=20, transform=ax[i].transAxes)
        
    ax[-1].set_xlabel("z [Mm]")
    plt.subplots_adjust(hspace=0.1)
    
    figname = "2D_Bhist.pdf"
    plt.savefig("figures/" + figname, bbox_inches="tight")
    
plot_2D_Bhist(5700, 500)