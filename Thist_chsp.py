import numpy as np
import matplotlib.pyplot as plt

#from scipy.interpolate import RegularGridInterpolator
#from matplotlib.ticker import FormatStrFormatter

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

HP_dict = np.load("files/HP_dict.pkl", allow_pickle=True)
TR_dict = np.load("files/TR_dict.pkl", allow_pickle=True)

Teff_list = [6500, 5770, 5000, 4000, 3200]
B0_list = [0,100,200,500]

ls_list = ["solid", "dashed", "dashdot", "dotted"]

fs = np.array([12,9])*0.8
fig, ax = plt.subplots(5,1, figsize=(fs[0],fs[1]), sharex=True, squeeze=True)

for i,Teff in enumerate(Teff_list):
    color = f"C{i}"
    label = r"$T_{eff}$ = "+f"{Teff} K"
    for j,B0 in enumerate(B0_list):
        F = h5py.File(input_files[Teff][B0], "r")

        logTR = TR_dict[Teff][B0]
        k_chsp = np.argmin(abs(logTR - (-5)))
        #print(k_chsp)

        bin_width = 25
        Trange = (1500, 7000)
        bins = int((Trange[-1] - Trange[0])/bin_width)

        T_hist, bin_edges = np.histogram(np.array(F["temperature"][0,:,:,:k_chsp]).flatten(), 
                                         range=Trange, bins=bins, density=True)
        
        if B0==0:
            ax[i].stairs(T_hist*100, bin_edges, color=color, ls=ls_list[j], lw=1.2, fill=False, label=label)
        else:
            ax[i].stairs(T_hist*100, bin_edges, color=color, ls=ls_list[j], lw=1.2, fill=False)
               
    
    #tls = ax[i].get_yticklabels()
    #tls[0] = ""
    #ax[i].set_yticklabels(tls)
    yticks = ax[i].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    
    ax[i].legend()
    ax[i].set_axisbelow(True)
    ax[i].grid()
    #ax[i].text(0.82, 0.7, r"$T_{eff}$" + f" = {Teff} K", transform=ax[i].transAxes)

fig.text(0.06, 0.5, "% of grid cells", ha="center", va="center", rotation="vertical")
ax[-1].set_xlabel("Gas temperature [K]")
plt.subplots_adjust(hspace=0)

figname = "T_hist_chsp.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")