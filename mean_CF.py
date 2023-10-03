import numpy as np
import matplotlib.pyplot as plt

import h5py

path_Lf3D = "/mn/stornext/d19/RoCS/alma/emissa_sim/linfor3D/outhdf/"
f = {"500nm" : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_05000_mu1_00_linfor_3D_2.hdf", "r"),
     "1mm"   : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_01mm_mu1_00_linfor_3D_2.hdf", "r"),
     "3mm"   : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_03mm_mu1_00_linfor_3D_2.hdf", "r")} 

F = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/tst/d3t57g44_v000G_n019_art_it000_mode1.h5", "r")

fig,ax = plt.subplots(figsize=(8,6))

for i,key in enumerate(f.keys()):
    
    z_CF = np.array(f[key]["contfctz"])*1e-8
    k = np.argmax(z_CF)
    z_CF = z_CF[:k+1]
    
    CF = np.mean(np.array(f[key]["contfunc"][:k+1,:,:]), axis=(1,2))
    
    ax.plot(z_CF, CF/np.max(CF))
    ax.fill_between(z_CF, np.zeros(len(CF)), CF/np.max(CF), alpha=0.5, label=r"$\lambda = $ {:} {:}".format(key[:-2], key[-2:]))
    
    ax.set_axisbelow(True)
    ax.grid()
    
    #ax.axvline(x=

ax.legend()
ax.set_xlabel("z [Mm]")
ax.set_ylabel("norm. cont. func.")

figname="mean_CF.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")