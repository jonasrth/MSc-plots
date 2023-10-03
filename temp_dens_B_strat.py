import numpy as np
import matplotlib.pyplot as plt

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

def plot_temp_dens(Teff):
    """
    """
    
    fig,ax = plt.subplots(2, figsize=(10,6), sharex=True)
    
    ### Adjust z origin to Rosseland mean opacity unity
    F = h5py.File(input_files[Teff][0], "r")
    f = h5py.File(output_files[Teff][0][0.0005], "r")
    
    z = np.array(F["z"][0,0,0,:])*1e-8
    z_top = np.max(np.array(f["contfctz"]))*1e-8
    z = z - z[0] + z_top
    
    sorted_keys = np.sort(list(input_files[Teff].keys()))
    for i,B in enumerate(sorted_keys):
        F = h5py.File(input_files[Teff][B], "r")
        z_ = np.array(F["z"][0,0,0,:])*1e-8
        z_ = z_ - z_[-1] + z[-1]
        #f = h5py.File(output_files[Teff][B][0.0005], "r")
        
        T = np.mean(np.array(F["temperature"][0]), axis=(0,1))
        rho = np.mean(np.array(F["dens"][0]), axis=(0,1))
        
        ax[0].plot(z_, T*1e-3, label=f"{B} G")
        ax[1].plot(z_, rho)
        
    handles, labels = ax[0].get_legend_handles_labels()
    order = - (1 + np.arange(len(handles)))
    ax[0].legend([handles[idx] for idx in order],[labels[idx] for idx in order])
    
    ax[0].grid()
    ax[0].set_ylabel("Gas temperature [kK]")
    
    ax[1].grid()
    ax[1].set_yscale("log")
    ax[1].set_ylabel("Mass density [g/cm$^3$]")
    
    ax[-1].set_xlabel("z [Mm]")
    
    plt.subplots_adjust(hspace=0.1)
    
    figname = "T_rho_B_strat.pdf"
    plt.savefig("figures/" + figname, bbox_inches="tight")
    
plot_temp_dens(5700)