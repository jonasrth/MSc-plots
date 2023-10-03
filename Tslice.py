import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.signal import find_peaks

import h5py

path = "/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/"

F = {
     6500 : h5py.File(path + "d3t65g45c_v000G_end_art.h5", "r"),
     5700 : h5py.File(path + "tst/d3t57g44_v000G_n019_art_it000_mode1.h5", "r"),
     5000 : h5py.File(path + "d3t50g45c_v000G_n009_end_art.h5", "r"),
     4000 : h5py.File(path + "d3t40g45c_v000G_n012_end_art.h5", "r"),
     3200 : h5py.File(path + "d3t32g45c_v000G_n004_it000_full_art.h5", "r")
    }

def surface_height_placeholder(atmosfile, f=0.4):

    Tstrat = np.mean(np.array(atmosfile["temperature"][0,...]), axis=(0,1))
    std = np.std(np.array(atmosfile["temperature"][0,...]), axis=(0,1))
    z = np.array(atmosfile["z"][0,0,0,:])

    rel_std = std/Tstrat

    k_phsp = find_peaks(rel_std, width=10)[0][-1]
    #k_Tmin = find_peaks(-rel_std, width=10)[0][-1]
    
    return z[k_phsp]

letters = ["a", "b", "c", "d", "e", "f"]

fig, ax = plt.subplots(3,2, figsize=(10,12))

"""
im = ax.flat[0].imshow(np.log10(T_new[:,::-1].T), 
                  extent=[0,np.sum(x)*1e-8,z_new[-1]*1e-3,z_new[0]*1e-3], vmin=np.log10(1500), vmax=np.log10(27000),
                  origin="lower", aspect="auto", cmap="inferno")
ax.flat[0].text(0.8, 0.8, letters[0], color="white", size=20, transform=ax.flat[0].transAxes)
"""

iy = 50
for i, Teff in enumerate(F.keys()):
        
    atmosfile = F[Teff]
    T = np.array(atmosfile["temperature"][0,:,iy,:])
    z = np.array(atmosfile["z"][0,0,0,:])
    dx = np.array(atmosfile["dx"])*1e-8
    
    if Teff==3200:
        # reshape T
        
        dz = 4*1e5
        Nz = int(abs(z[-1] - z[0])/dz)

        z_top = np.cumsum(np.ones(Nz + 1)*dz) - dz
        z_new = z[0] - z_top

        T = np.array(F[Teff]["temperature"][0,0,:,:])
        T_new = np.zeros((T.shape[0], len(z_new)))

        for j in range(T.shape[0]):
            f = interp1d(z, T[j,:])
            T_new[j,:] = f(z_new)
            
        T = T_new
        z = z_new

    z_s = surface_height_placeholder(atmosfile)
    z = z - z_s

    im = ax.flat[i].imshow(np.log10(T[:,::-1].T), 
                          extent=[0, np.sum(dx), z[-1]*1e-8, z[0]*1e-8], vmin=np.log10(1500), vmax=np.log10(27000),
                          aspect="auto", origin="lower", cmap="inferno")
    ax.flat[i].text(0.8, 0.8, letters[i], color="white", size=20, transform=ax.flat[i].transAxes)

for j in range(3):
    ax[j,0].set_ylabel("z [Mm]")
for i in range(2):
    ax[-1,i].set_xlabel("x [Mm]")
    
plt.subplots_adjust(hspace=0.1)
fig.colorbar(im, ax=ax.ravel().tolist(), location="top", fraction=0.046, pad=0.02, label=r"$\log(T$ [K]$)$")

figname = "Tslice.pdf"
plt.savefig("figures/" + figname, bbox_inches="tight")