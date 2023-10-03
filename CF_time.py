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

def plot_CF_time(inds, t_list, figname="CF_time.pdf"):
    """
    """
    ix,iy = inds[0],inds[1]
    
    zlim = [0.3,1.4]
    
    z = np.array(F[0]["z"][0,0,0,:])*u.cm
    
    zinds = [np.argmin(abs(z.to(u.Mm).value - zlim[1])), np.argmin(abs(z.to(u.Mm).value - zlim[0]))+1]
    
    z = z[zinds[0]:zinds[1]]
    
    #CF = np.array(f[5]["CFunc"][0,ix,iy,:])
    #CF = CF/np.max(CF)
    
    fig,ax = plt.subplots(3, figsize=(10,9), sharex=True)
    
    ls_list = ["solid", "dashed", "dashdot", "dotted"]
    
    for i,n in enumerate(t_list):
        
        CF = np.array(f[n]["CFunc"][0,ix,iy,zinds[0]:zinds[1]])
        CF = CF/np.max(CF)
        
        T = np.array(F[n]["temperature"][0,ix,iy,zinds[0]:zinds[1]])
        
        ne = np.array(F[n]["xne"][0,ix,iy,zinds[0]:zinds[1]])
        
        label = f"t = {time_list[n]} s"
        ax[0].plot(z.to(u.Mm), CF, color="black", ls=ls_list[i%len(ls_list)], label=label)
        
        ax[1].plot(z.to(u.Mm), T, color="black", ls=ls_list[i%len(ls_list)])
        
        ax[2].plot(z.to(u.Mm), ne, color="black", ls=ls_list[i%len(ls_list)])
    
    for i in range(len(ax)):
        ax[i].grid()
        
    ax[0].legend()
    ax[0].set_ylabel("CF (normalised)")
    
    ax[1].set_ylabel("Gas temperature [K]")
    
    ax[2].set_ylabel(r"Electron density [cm$^{-3}$]")
    ax[2].set_yscale("log")
    
    ax[-1].set_xlim(zlim)
    ax[-1].set_xlabel("z [Mm]")
    
    plt.savefig("figures/"+figname, bbox_inches="tight")
    
ix,iy = 26, 225
inds = [ix,iy]

plot_CF_time(inds, [1, 5, 9, 17], figname="CF_time.pdf")