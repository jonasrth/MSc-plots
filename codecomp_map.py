import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import cumtrapz

from astropy import units as u
from astropy import constants as const

I_units_rh = u.W*u.m**(-2)*u.Hz**(-1)*u.sr**(-1)
I_units_art = u.erg*u.cm**(-2)*u.s**(-1)*u.Hz**(-1)*u.sr**(-1)
I_units_Lf = u.erg*u.cm**(-2)*u.s**(-1)*u.angstrom**(-1)*u.sr**(-1)

c = const.c
h = const.h
kB = const.k_B

import h5py

import sys

sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

#Lf3D_wave_list = ["5000", "00_30", "01_00", "03_00", "06_00", "09_00"]
art_wave_list = ["500nm", "00_30mm", "01_00mm", "03_00mm", "06_00mm", "09_00mm"]
wave_list = [0.0005, 0.3, 1, 3, 6, 9]

art_path = "/mn/stornext/d19/RoCS/jonasrth/ART/code_comparison/"
rh_path = "/mn/stornext/d19/RoCS/jonasrth/RH/code_comparison/d3t57g44c_v000G_n013_it000/"


F = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/d3t57g44c_v000G_n013_it000_full_art.h5", "r")
z = np.array(F["z"][0,0,0,:])*u.cm

def brightness_temperature(intensity,wave):
    """
    intensity in cgs units (per frequency)
    wavelength in Ångström
    """
    Tb = (h*c/(wave*kB)*np.log(1 + (2*h*c/(wave**3*intensity*u.sr)).to(1))**(-1)).to("K")
    return Tb



wave = 3*u.mm

fs = np.array([15,10])*0.6
fig,ax = plt.subplots(2,3,figsize=(fs[0],fs[1]))

DX = np.sum(np.array(F["dx"]))*1e-8
extent = [0,DX,0,DX]

f_L = h5py.File(output_files[5770][0][3], "r")
I_L = np.array(f_L["intensity"])*I_units_Lf
Tb_L = brightness_temperature((wave**2/c)*I_L, wave).to(u.kK)

i = 3
f_A = h5py.File(art_path + "d3t57g44_000G_CF_" + art_wave_list[i] + ".h5", "r")
I_A = np.array(f_A["Stokes_I"][0,:,:,0])*I_units_art
Tb_A = brightness_temperature(I_A, wave).to(u.kK)

f_R = h5py.File(rh_path + "output_ray.hdf5", "r")
I_R = np.array(f_R["intensity"][...,-3])*I_units_rh
Tb_R = brightness_temperature(I_R, wave).to(u.kK)

vmin = np.min([Tb_L, Tb_A, Tb_R])
vmax = np.max([Tb_L, Tb_A, Tb_R])

ax[0,0].imshow(Tb_L, origin="lower", extent=extent, vmin=vmin, vmax=vmax, cmap="inferno")
ax[0,0].set_title("Linfor3D")
ax[0,0].set_ylabel("y [Mm]")

ax[0,1].imshow(Tb_A, origin="lower", extent=extent, vmin=vmin, vmax=vmax, cmap="inferno")
ax[0,1].set_title("ART")

im0 = ax[0,2].imshow(Tb_R, origin="lower", extent=extent, vmin=vmin, vmax=vmax, cmap="inferno")
ax[0,2].set_title("RH")

d_LR = Tb_L - Tb_R
d_LA = Tb_L - Tb_A

vmin = np.min([d_LA, d_LR])
vmax = np.max([d_LA, d_LR])

ax[1,0].imshow(d_LR, origin="lower", extent=extent,  vmin=vmin, vmax=vmax, cmap="rainbow")
ax[1,0].set_title("Linfor3D - RH")
ax[1,0].set_xlabel("x [Mm]")
ax[1,0].set_ylabel("y [Mm]")

im1 = ax[1,1].imshow(d_LA, origin="lower", extent=extent, vmin=vmin, vmax=vmax, cmap="rainbow")
ax[1,1].set_title("Linfor3D - ART")
ax[1,1].set_xlabel("x [Mm]")

plt.subplots_adjust(wspace=0.25)

fig.colorbar(im0, ax=ax[0,:].ravel().tolist(), location="left", fraction=0.046, pad=0.06, label="Brightness temperature [kK]")
fig.colorbar(im1, ax=ax[1,:].ravel().tolist(), location="left", fraction=0.046, pad=0.06, label="Absolute error [kK]")


H_L, L_edges = np.histogram(Tb_L.value.flatten(), density=True, bins=100)
H_A, A_edges = np.histogram(Tb_A.value.flatten(), density=True, bins=100)
H_R, R_edges = np.histogram(Tb_R.value.flatten(), density=True, bins=100)

ax[1,2].stairs(H_L, L_edges, label="Linfor3D")
ax[1,2].stairs(H_A, A_edges, label="ART")
ax[1,2].stairs(H_R, R_edges, label="RH")

ax[1,2].legend()
ax[1,2].set_title("Histogram")
ax[1,2].yaxis.tick_right()
ax[1,2].yaxis.set_label_position("right")
ax[1,2].set_title("")
ax[1,2].set_xlabel("Brightness temperature [kK]")
ax[1,2].set_ylabel("Occurence")

letter_list = ["a","b","c","d","e","f"]
for i in range(5):
    ax.flat[i].text(0.1,0.9, letter_list[i], size = 15, color="white", transform=ax.flat[i].transAxes)
    ax.flat[i].set_xlabel("x [Mm]")
    ax.flat[i].set_ylabel("y [Mm]")
ax.flat[-1].text(0.1,0.9, letter_list[-1], size = 15, color="black", transform=ax.flat[-1].transAxes)

#plt.tight_layout()
figname="codecomp_maps.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")