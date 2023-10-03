import numpy as np
import matplotlib.pyplot as plt

import h5py

Linfor3D_path = "/mn/stornext/d19/RoCS/alma/emissa_sim/linfor3D/outhdf/"

"""
# 500 nm
f = {}
f[0] = h5py.File(Linfor3D_path + "d3t57g44c_v000G_n019_it000_05000_mu1_00_linfor_3D_2.hdf", "r")
f[100] = h5py.File(Linfor3D_path + "d3gt57g44v100_chro01_n016_it003_05000_mu1_00_linfor_3D_2.hdf", "r")
f[200] = h5py.File(Linfor3D_path + "d3gt57g44v200_chro01_n016_it001_05000_mu1_00_linfor_3D_2.hdf", "r")
f[500] = h5py.File(Linfor3D_path + "d3gt57g44v500_chro01_n016_it001_05000_mu1_00_linfor_3D_2.hdf", "r")
"""

# 3 mm
f = {}
f[0] = h5py.File(Linfor3D_path + "d3t57g44c_v000G_n019_it000_03mm_mu1_00_linfor_3D_2.hdf", "r")
f[100] = h5py.File(Linfor3D_path + "d3gt57g44v100_chro01_n016_it003_03_00mm_mu1_00_linfor_3D_2.hdf", "r")
f[200] = h5py.File(Linfor3D_path + "d3gt57g44v200_chro01_n016_it001_03_00mm_mu1_00_linfor_3D_2.hdf", "r")
f[500] = h5py.File(Linfor3D_path + "d3gt57g44v500_chro01_n016_it001_03_00mm_mu1_00_linfor_3D_2.hdf", "r")


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

fig,ax = plt.subplots(figsize=(10,6), sharex=True)

bin_width = 0.01
zrange = (0.2, 1.5)
bins = int((zrange[-1] - zrange[0])/bin_width)

for i,G in enumerate(f.keys()):
    
    zt1 = find_mean_z_CF(f[G])
    
    H,bin_edges = np.histogram(zt1.flatten(), bins=bins, range=zrange)
    
    ax.stairs(H/np.max(H), bin_edges, color=f"C{i}", fill=True, label=fr"$B_0$ = {G} G", alpha=0.7)
    ax.stairs(H/np.max(H), bin_edges, color=f"C{i}", lw=2)

ax.text(0.065, 0.9, "$\lambda = 3.0$ mm", size=12, transform=ax.transAxes)
    
ax.legend()
    
ax.set_axisbelow(True)
ax.grid()
ax.set_xlabel("z [Mm]")
ax.set_ylabel("Occurence (normalised)")

figname = "zt1_hist_mag.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")