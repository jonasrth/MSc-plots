import numpy as np
import matplotlib.pyplot as plt

import h5py

from astropy import constants as const
from astropy import units as u

c = const.c
h = const.h
kB = const.k_B

I_units_Lf3D = u.erg*u.cm**(-2)*u.s**(-1)*u.angstrom**(-1)*u.sr**(-1)
I_units = u.erg*u.cm**(-2)*u.s**(-1)*u.Hz**(-1)*u.sr**(-1)


path_Lf3D = "/mn/stornext/d19/RoCS/alma/emissa_sim/linfor3D/outhdf/"
f = {"500nm" : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_05000_mu1_00_linfor_3D_2.hdf", "r"),
          "1mm"   : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_01mm_mu1_00_linfor_3D_2.hdf", "r"),
          "3mm"   : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_03mm_mu1_00_linfor_3D_2.hdf", "r")} 

F = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/tst/d3t57g44_v000G_n019_art_it000_mode1.h5", "r")

def brightness_temperature(intensity,wave):
    """
    intensity in cgs units (per frequency)
    wavelength in Ångström
    """
    #Tb = (wave**2/(2*kB)*intensity*u.sr).to("K")
    Tb = (h*c/(wave*kB)*np.log(1 + (2*h*c/(wave**3*intensity*u.sr)).to(1))**(-1)).to("K")
    return Tb

dx,dy = np.array(F["dx"][0])*1e-8, np.array(F["dy"][0])*1e-8

shape = f["500nm"]["intensity"].shape
extent = [0 - dx/2, shape[0]*dx + dx/2, 0 - dy/2, shape[1]*dy + dy/2]

fig,ax = plt.subplots(1,3, figsize=(15,5), sharey=True)

text_list = ["a", "b", "c"]
for i,key in enumerate(f.keys()):
    wave = np.array(f[key]["wavelength"])*u.nm
    intensity = ((wave**2/c)*np.array(f[key]["intensity"])*I_units_Lf3D).to(I_units)
    #intensity = np.array(f[key]["intensity"])
    im = ax[i].imshow(brightness_temperature(intensity, wave), origin="lower", extent=extent)
    plt.colorbar(im, ax=ax[i], location="top", fraction=0.046, pad=0.04, label="Gas temperature [K]")
    
    ax[i].text(0.9, 0.9, text_list[i], size=20, color="white", transform=ax[i].transAxes)
    
    ax[i].set_xlabel("x [Mm]")

ax[0].set_ylabel("y [Mm]")
plt.subplots_adjust(wspace=0)

figname = "int_img.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")