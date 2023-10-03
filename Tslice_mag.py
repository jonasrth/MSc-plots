import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import find_peaks

import h5py

input_path = "/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/"

F = {}
F["000G"] = h5py.File(input_path + "d3t50g45c_v000G_n009_end_art.h5", "r")
F["100G"] = F["000G"]
F["200G"] = h5py.File(input_path + "d3t50g45c_v200G_n002_it000_full_art.h5", "r")
F["500G"] = h5py.File(input_path + "d3t50g45c_v500G_n002_it000_full_art.h5", "r")

def surface_height_placeholder(atmosfile, f=0.4):

    Tstrat = np.mean(np.array(atmosfile["temperature"][0,...]), axis=(0,1))
    std = np.std(np.array(atmosfile["temperature"][0,...]), axis=(0,1))
    z = np.array(atmosfile["z"][0,0,0,:])

    rel_std = std/Tstrat

    k_phsp = find_peaks(rel_std, width=10)[0][-1]
    #k_Tmin = find_peaks(-rel_std, width=10)[0][-1]
    
    return z[k_phsp]

z = np.array(F["000G"]["z"][0,0,0,:])*1e-8
L = np.cumsum(np.array(F["000G"]["dx"]))*1e-8

iy = 50
figname = "Tslice_mag_t50.pdf"

fig,ax = plt.subplots(2,2, figsize=(10,7), sharey=True, sharex=True)

#T = {}
txt_list = ["a", "b", "c", "d"]
for i,G in enumerate(F.keys()):
    T = np.array(F[G]["temperature"][0,:,iy,:])[:,::-1].T
    z_new = z - surface_height_placeholder(F[G])*1e-8
    im = ax.flat[i].imshow(np.log10(T), origin="lower", extent=[L[0], L[-1], z_new[-1], z_new[0]], cmap="inferno")
    ax.flat[i].text(0.8, 0.8, txt_list[i], color="white", size=20, transform=ax.flat[i].transAxes)

ax[0,0].set_ylabel("z [Mm]")
ax[1,0].set_ylabel("z [Mm]")

ax[1,0].set_xlabel("y [Mm]")
ax[1,1].set_xlabel("y [Mm]")

plt.subplots_adjust(wspace=0, hspace=0)

fig.colorbar(im, ax=ax.ravel().tolist(), fraction=0.046, pad=0.04, label=r"$\log(T$ [K]$)$")

plt.savefig("figures/" + figname, bbox_inches="tight")