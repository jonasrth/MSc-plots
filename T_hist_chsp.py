import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import find_peaks

import h5py

path = "/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/"

F = h5py.File(path + "tst/d3t57g44_v000G_n019_art_it000_mode1.h5", "r")
F = h5py.File(path + "d3t50g45c_v000G_n009_end_art.h5", "r")
z = np.array(F["z"][0,0,0,:])

F = {
     6500 : h5py.File(path + "d3t65g45c_v000G_end_art.h5", "r"),
     5700 : h5py.File(path + "tst/d3t57g44_v000G_n019_art_it000_mode1.h5", "r"),
     5000 : h5py.File(path + "d3t50g45c_v000G_n009_end_art.h5", "r"),
     4000 : h5py.File(path + "d3t40g45c_v000G_n012_end_art.h5", "r"),
     3200 : h5py.File(path + "d3t32g45c_v000G_n004_it000_full_art.h5", "r")
    }

def find_chsp_index(atmosfile, f=0.5):

    Tstrat = np.mean(np.array(atmosfile["temperature"][0,...]), axis=(0,1))
    std = np.std(np.array(atmosfile["temperature"][0,...]), axis=(0,1))

    rel_std = std/Tstrat

    k_phsp = find_peaks(rel_std, width=10)[0][-1]
    k_Tmin = find_peaks(-rel_std, width=10)[0][-1]
    
    rel_std_diff = abs(rel_std[k_phsp] - rel_std[k_Tmin])

    k_chsp = np.argwhere(rel_std < rel_std[k_Tmin] + f*rel_std_diff)[0][0]
    
    return k_chsp

fig, ax = plt.subplots(5,1, figsize=(12,9), sharex=True, sharey=True)


for i,Teff in enumerate(F.keys()):
    atmosfile = F[Teff]
    
    k_chsp = find_chsp_index(atmosfile, f = 0.4)
    
    bin_width = 25
    Trange = (1500, 7000)
    bins = int((Trange[-1] - Trange[0])/bin_width)
    
    T_hist, bin_edges = np.histogram(np.array(atmosfile["temperature"][0,:,:,:k_chsp]).flatten(), 
                                     range=Trange, bins=bins)
    
    ax[i].set_axisbelow(True)
    ax[i].grid()
    ax[i].stairs(T_hist/np.max(T_hist), bin_edges, fill=True)
    ax[i].text(0.82, 0.7, r"$T_{eff}$" + f" = {Teff} K", transform=ax[i].transAxes)

fig.text(0.06, 0.5, "Normalised occurence", ha="center", va="center", rotation="vertical")
ax[-1].set_xlabel("Gas temperature [K]")
plt.subplots_adjust(hspace=0.1)

figname = "T_hist_chsp.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")