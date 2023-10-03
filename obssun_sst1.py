import numpy as np
import matplotlib.pyplot as plt

import h5py

F = h5py.File("/mn/stornext/d19/RoCS/svenwe/lecture/AST5770/data/obssun_sst1/obssun_sst1.h5", "r")

I = np.mean(np.array(F["intensity"][0,...]), axis=0)

x = np.array(F["x"])
y = np.array(F["y"])

fig,ax = plt.subplots(figsize=(5,5))

extent = [x[0],x[-1],y[0],y[-1]]
im = ax.imshow(I, origin="lower", extent=extent, cmap="inferno")
fig.colorbar(im, fraction=0.047, pad=0.02, label="Intensity [counts]")

ax.set_xlabel("X [arcsec]")
ax.set_ylabel("Y [arcsec]")

figname="obssun_sst1.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")