import numpy as np
import matplotlib.pyplot as plt

import h5py

from astropy import constants
from astropy import units as u

import h5py

kB = constants.k_B
c = constants.c
m_p = constants.m_p
m_e = constants.m_e

u_int = u.erg*u.cm**(-2)*u.s**(-1)*u.Hz**(-1)

input_path = "/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/"
output_path = "/mn/stornext/d19/RoCS/jonasrth/ART/time_series_CF/"
input_path = "/mn/stornext/d19/RoCS/jonasrth/ART_input/time_series_backup/"

n_list = [13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19]
it_list = [0, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]
time_list = [1800.01, 1810.04, 1820.00, 1830.02, 1840.00, 1850.02, 1860.00, 1870.01, 1880.01, 1890.02, 
             1900.00, 1910.01, 1920.02, 1930.00, 1940.01, 1950.01, 1960.03, 1970.01, 1980.02]

f = {}
F = {}
for i in range(len(n_list)):
    f[i] = h5py.File(output_path + f"output_n0{n_list[i]}_it00{it_list[i]}.h5", "r")
    F[i] = h5py.File(input_path + f"d3t57g44c_v000G_n0{n_list[i]}_it00{it_list[i]}_full_art.h5", "r")

wave = f[0]["Wavelength"][0]*u.angstrom



def plot_CF_time_image(inds, figname="CF_time_image.pdf"):
    """
    """
    ix,iy = inds[0], inds[1]
    
    zlim = [0.3,1.4]
    
    z = np.array(F[0]["z"][0,0,0,:])*u.cm
    
    zinds = [np.argmin(abs(z.to(u.Mm).value - zlim[1])), np.argmin(abs(z.to(u.Mm).value - zlim[0]))+1]
    
    z = z[zinds[0]:zinds[1]]
    
    CF_im = np.zeros((len(z), len(n_list)))
    zt1 = np.zeros(len(n_list))
    for i,n in enumerate(n_list):
        CF = np.array(f[i]["CFunc"][0,ix,iy,zinds[0]:zinds[1]])
        CF_im[:,i] = CF[::-1]
        zt1[i] = np.array(f[i]["Tau1"][0,ix,iy,0])*1e-8
        
    
    fig,ax = plt.subplots(figsize=(10,6))
    
    dt = abs(time_list[1]-time_list[0])
    dz = abs(z[0] - z[1]).to(u.Mm).value
    extent = [time_list[0] - dt/2, time_list[-1] + dt/2, z[-1].to(u.Mm).value - dz/2, z[0].to(u.Mm).value + dz/2]
    
    ax.set_axisbelow(True)
    ax.grid(alpha=0.5)
    
    im = ax.imshow(CF_im, cmap="Greys", aspect="auto", interpolation="nearest", extent=extent, origin="lower")
    plt.colorbar(im, label=r"Contribution function [10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ cm$^{-1}$]")
    
    ax.plot(time_list, zt1, ls="", marker="o", label=r"$z(\tau=1)$")
    
    ax.legend()
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("z [Mm]")
    
    plt.savefig("figures/" + figname, bbox_inches="tight")
    

def plot_time_quantity(inds, quantity="T", figname="T_time_image.pdf"):
    """
    """
    
    ix,iy = inds[0], inds[1]
    
    zlim = [0.3,1.4]
    
    z = np.array(F[0]["z"][0,0,0,:])*u.cm
    
    zinds = [np.argmin(abs(z.to(u.Mm).value - zlim[1])), np.argmin(abs(z.to(u.Mm).value - zlim[0]))+1]
    
    z = z[zinds[0]:zinds[1]]
    
    if quantity=="T":
        var   = "temperature"
        scale = lambda x: x
        label = r"Gas Temperature [K]"
    elif quantity=="ne":
        var   =  "xne"
        label = r"$\log_{10}($Electron number density [1/cm$^3$]$)$"
        scale = lambda x: np.log10(x)
        
    data = np.zeros((len(z), len(n_list)))
    zt1 = np.zeros(len(n_list))
    for i,n in enumerate(n_list):
        
        col = np.array(F[i][var][0,ix,iy,zinds[0]:zinds[1]])
        data[:,i] = scale(col[::-1])
        
        # Some columns of ne off with factor of 1000 for some reason
        if quantity=="ne":
            if i==10 or i==16:
                data[:,i] += -3
        
        zt1[i] = np.array(f[i]["Tau1"][0,ix,iy,0])*1e-8
        
    fig,ax = plt.subplots(figsize=(10,6))
    
    dt = abs(time_list[1]-time_list[0])
    dz = abs(z[0] - z[1]).to(u.Mm).value
    extent = [time_list[0] - dt/2, time_list[-1] + dt/2, z[-1].to(u.Mm).value - dz/2, z[0].to(u.Mm).value + dz/2]
    
    im = ax.imshow(data, cmap="inferno", aspect="auto", interpolation="nearest", extent=extent, origin="lower")
    plt.colorbar(im, label=label)
    
    ax.plot(time_list, zt1, ls="", marker="o", label=r"$z(\tau=1)$")
    
    ax.set_axisbelow(True)
    ax.grid(alpha=0.5)
    
    ax.legend()
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("z [Mm]")
    
    plt.savefig("figures/" + figname, bbox_inches="tight")
    
ix,iy = 26, 225
inds = [ix,iy]

plot_CF_time_image(inds, figname = "CF_time_image.pdf")
#plot_time_quantity(inds, quantity = "T",  figname = "T_time_image.pdf")
#plot_time_quantity(inds, quantity = "ne", figname = "ne_time_image.pdf")