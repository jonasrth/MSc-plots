import numpy as np
import matplotlib.pyplot as plt

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

"""
f = {
    "500nm" : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/test_mode1_500nm.h5", "r"),
    "1mm"   : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/test_mode1_1mm.h5", "r"),
    "3mm"   : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/test_mode1_3mm.h5", "r")
    }

f = {
    "500nm" : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/new_mode1_500nm.h5", "r"),
    "1mm"   : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/new_mode1_1mm.h5", "r"),
    "3mm"   : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/new_mode1_3mm.h5", "r")
    }

path_Lf3D = "/mn/stornext/d19/RoCS/alma/emissa_sim/linfor3D/outhdf/"
f_Lf3D = {"500nm" : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_05000_mu1_00_linfor_3D_2.hdf", "r"),
          "1mm"   : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_01mm_mu1_00_linfor_3D_2.hdf", "r"),
          "3mm"   : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_03mm_mu1_00_linfor_3D_2.hdf", "r")} 

#F = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/d3t57g44c_v000G_n019_it000_full_art.h5", "r")
F = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/tst/d3t57g44_v000G_n019_art_it000_mode1.h5", "r")
"""

def find_mean_z_CF(linfor3D_output):
    """
    """
    
    z_CF = np.array(linfor3D_output["contfctz"])*1e-8
    k = np.argmax(z_CF)
    z_CF = z_CF[:k+1]
    
    CF = np.array(linfor3D_output["contfunc"][:k+1,:,:])
    
    #z_tiled = np.tile(z_CF, (CF.shape[1], CF.shape[2], 1)).T
    #z = np.average(z_tiled, weights=CF, axis=0)
    
    mu = np.trapz(z_CF[:,np.newaxis,np.newaxis]*CF, z_CF, axis=0)\
    /np.trapz(CF, z_CF, axis=0)
    
    return mu

def find_height_of_formation(linfor3D_output, f=0.5):
    """
    Finds proxy for height of formation from contribution function (CF) given by linfor3D. Found by crude interpolation
    between points closest to height where fraction f is reached.
    
    f - [float] fraction f of CF comes from above output height
    """
    
    # Will be integrating from bottom up
    f = 1 - f
    
    # Restricting to height range where CF is recorded by Linfor3D
    z_CF = np.array(linfor3D_output["contfctz"])*1e-8
    k = np.argmax(z_CF)
    z_CF = z_CF[:k+1]
    
    CF = np.array(linfor3D_output["contfunc"][:k+1,:,:])
    
    CF_cumsum = np.cumsum(CF, axis=0)
    CF_cumsum_norm = CF_cumsum/CF_cumsum[-1]
    
    #interpolating:
    
    inds1 = np.argmin(abs(CF_cumsum_norm - f), axis=0)
    f1 = np.take_along_axis(CF_cumsum_norm, inds1[np.newaxis, ...], axis=0)[0,...]
    z1 = z_CF[inds1]
    
    inds2 = np.where(f > f1, inds1+1, inds1-1)
    f2 = np.take_along_axis(CF_cumsum_norm, inds2[np.newaxis, ...], axis=0)[0,...]
    z2 = z_CF[inds2]
    
    deltaz = abs(z1 - z2)
    deltaf = abs(f1 - f2)
    f_frac = f - f1
    
    z_star = z1 + (f_frac/deltaf)*deltaz
    
    return z_star

fig,ax = plt.subplots(figsize=(10,6), sharex=True)

bin_width = 0.01
zrange = (-0.25, 1.45)
bins = int((zrange[-1] - zrange[0])/bin_width)

Teff = 5700
for i,key in enumerate(output_files[Teff][0].keys()):
    
    #zt1 = np.array(f[key]["Tau1"][0,:,:,0])*1e-8
    #zt1 = find_height_of_formation(f_Lf3D[key])
    zt1 = find_mean_z_CF(h5py.File(output_files[Teff][0][key], "r"))
    
    zt1_mean = np.mean(zt1).astype(np.float64)
    zt1_std = np.std(zt1).astype(np.float64)
    
    H,bin_edges = np.histogram(zt1.flatten(), bins=bins, range=zrange)
    
    ax.stairs(H/np.max(H), bin_edges, color=f"C{i}", fill=True, label=r"$\lambda = $ {:} mm".format(key), alpha=0.7)
    ax.stairs(H/np.max(H), bin_edges, color=f"C{i}", lw=2)
    #ax[i].axvspan(zt1_mean - zt1_std, zt1_mean + zt1_std, color="grey", alpha=0.2)
    #ax[i].axvline(x=zt1_mean, ls="dashed", color="red", label=r"$\mu = {:.2f}$ Mm".format(zt1_mean))
    #ax[i].axvline(x=zt1_mean-zt1_std, ls="dashed", color="grey")
    #ax[i].axvline(x=zt1_mean+zt1_std, ls="dashed", color="grey", label=r"$\sigma = {:.2f}$ Mm".format(zt1_std))
    
    #ax[i].text(0.8, 0.8, fr"$\lambda = $ {key[:-2]} {key[-2:]}", transform=ax[i].transAxes)

"""
ax[0].set_ylabel("occurence (norm)")
ax[0].legend(loc="lower right")

z_CF = np.array(f_Lf3D["3mm"]["contfctz"])*1e-8
k = np.argmax(z_CF)
z_CF = z_CF[:k+1]

T = np.mean(np.array(F["temperature"][0,...]),axis=(0,1))
rho = np.mean(np.array(F["dens"][0,...]),axis=(0,1))

ax[1].plot(z_CF[::-1], T[:k+1])
ax[1].set_ylabel("Gas temperature [K]")

ax[2].plot(z_CF[::-1], rho[:k+1])
ax[2].set_yscale("log")
ax[2].set_ylabel(r"Mass density [g/cm$^3$]")
"""
ax.legend()
ax.set_axisbelow(True)
ax.grid()
ax.set_xlabel("z [Mm]")
ax.set_ylabel("Normalised occurence")
plt.subplots_adjust(hspace=0.1)

figname = "zt1_hist.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")