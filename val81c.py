import numpy as np
import matplotlib.pyplot as plt

import h5py

F = h5py.File("/mn/stornext/d19/RoCS/svenwe/lecture/AST5770/data/ref/val81c.h5", "r")

h = np.array(F["h"])*1e-3
rho = np.array(F["rho"])
T = np.array(F["t"])*1e-3

fig,ax = plt.subplots(figsize=(6,5))

ax.plot(h, rho, color="k", ls="dashed")

ax.set_yscale("log")
ax.set_xlabel("Height [Mm]")
#ax.grid(which="both")
ax.set_ylabel(r"Mass density [g/cm$^3$]")

ax2 = ax.twinx()
ax2.plot(h, T, color="k")
ax2.set_yscale("log")
ax2.set_ylabel("Gas temperature [kK]")

figname="val81c.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")