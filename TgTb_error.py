import matplotlib.pyplot as plt
import numpy as np

import pickle

test_dict = np.load("files/DT_dict.pkl", allow_pickle=True)

fs = np.array([10,7])*0.8
fig,ax = plt.subplots(figsize=(fs[0],fs[1]))

wave = test_dict["wave"]
ls_list = ["solid", "dashed", "dashdot", "dotted"]
for i,Teff in enumerate([6500,5770,5000,4000,3200]):
    color=f"C{i}"
    for j,B0 in enumerate([0,100,200,500]):
        
        dT_abs = test_dict[Teff][B0]["dT_abs"]
        dT_rel = test_dict[Teff][B0]["dT_rel"]
        
        if B0==0:
            ax.plot(wave,dT_rel, color=color, ls=ls_list[j], label=r"$T_{eff}$ = " + f"{Teff} K")
        else:
            ax.plot(wave,dT_rel, color=color, ls=ls_list[j])
ax.legend()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Wavelength [mm]")
ax.set_ylabel(r"$\langle|T_B - T_G(\tau=1)|/T_B\rangle$")
ax.set_xlim([wave[0], wave[-1]])
#ax.grid(which="both")

figname="Tg_Tb_error.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")



