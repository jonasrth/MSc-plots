import numpy as np
import matplotlib.pyplot as plt

import csv
import h5py

from astropy import constants as const
from astropy import units as u

data_path = "/mn/stornext/d20/RoCS/atulm/Project1_stars/SED_data/"

logg_dict = {"A": 4.31, "B": 4.20, "C": 4.35, "D": 4.25, "E": 4.25, "F": 4.30,
             "G": 4.29, "H": 4.42, "I": 4.53, "J": 4.60, "K": 4.69, "L": 5.16,
             "S": 4.44}

with open(data_path + "star_props.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    csv_array = np.array(list(csv_reader))
    
    #print(csv_array[:,-1])
    
    n_ID = np.argwhere(csv_array[0] == "ID")[0][0]
    n_star = np.argwhere(csv_array[0] == "Star")[0][0]
    n_Teff = np.argwhere(csv_array[0] == " T (K)(Catalogs)")[0][0]
    n_SpT = np.argwhere(csv_array[0] == "Sp. Type")[0][0]
    
    IDs = csv_array[1:, n_ID]
    stars = csv_array[1:, n_star]
    Teff = csv_array[1:, n_Teff].astype(np.float64)
    SpT = csv_array[1:, n_SpT]
    
    logg = np.array([logg_dict[IDs[i]] for i in range(len(IDs))])

logg_mod = np.array([4.5, 4.44, 4.5, 4.5, 4.5])
Teff_mod = np.array([6500, 5770, 5000, 4000, 3200])

fig,ax = plt.subplots(figsize=(9,7))
ax.plot(logg, Teff, "bo", ms=10)
ax.plot(logg_mod, Teff_mod, "ro", ms=10)
ax.axvline(4.44, color="grey", alpha=0.5, ls="dashed")
ax.axhline(5770, color="grey", alpha=0.5, ls="dashed")

ax.legend(["obs", "mod"], loc="upper center", ncol=2)
ax.set_xlabel("$\log_{10}($ Surface gravity [cgs])")
ax.set_ylabel("Effective Temperature [K]")


IDs_sorted = np.sort(IDs)
inds = np.argsort(IDs)
stars_sorted = stars[inds]
SpT_sorted = SpT[inds]

for i,star in enumerate(stars_sorted):
    star_name = star
    star_name = star_name.replace("Alf", r"$\alpha$")
    star_name = star_name.replace("Gam", r"$\gamma$")
    star_name = star_name.replace("Eta", r"$\eta$")
    star_name = star_name.replace("Eps", r"$\epsilon$")
    stars_sorted[i] = star_name

for i, txt in enumerate(IDs):
    ax.annotate(txt, (logg[i]+0.01, Teff[i]+100))
    ax.text(0.75, 0.95 - 0.05*i, f"{IDs_sorted[i]} : {stars_sorted[i]} ({SpT_sorted[i]})", transform=ax.transAxes)
    
figname = "plot_params.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")