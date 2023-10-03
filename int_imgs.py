import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from astropy import constants as const
from astropy import units as u

kB = const.k_B
c = const.c
h = const.h

I_unitsw = u.erg*u.cm**(-2)*u.s**(-1)*u.angstrom**(-1)*u.sr**(-1)

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files


ls_list = ["solid", "dashed", "dashdot", "dotted"]

def brightness_temperature(intensity,wave):
    """
    intensity in cgs units (per frequency)
    wavelength in Ångström
    """
    Tb = (h*c/(wave*kB)*np.log(1 + (2*h*c/(wave**3*intensity*u.sr)).to(1))**(-1)).to("K")
    return Tb

def plot_int_imgs(wave, figname, label, dT=0.35):

    fs = np.array([12,15])*0.5
    fig,ax = plt.subplots(5,4,figsize=(fs[0],fs[1]))

    Teff_list = [6500,5770,5000,4000,3200]
    B0_list = [0,100,200,500]

    for i,Teff in enumerate(Teff_list):
        for j, B0 in enumerate(B0_list):
            f = h5py.File(output_files[Teff][B0][wave], "r")
            F = h5py.File(input_files[Teff][B0], "r")
            x = np.cumsum(np.array(F["dx"]))*1e-8

            intensity = ((wave*u.mm)**2/c)*np.array(f["intensity"])*I_unitsw
            
            #Tb = brightness_temperature(intensity, wave*u.mm)

            I = intensity/np.mean(intensity)

            ax[i,j].get_xaxis().set_ticks([])
            if j!=0:
                ax[i,j].get_yaxis().set_ticks([])
            
            extent = [x[0],x[-1],x[0],x[-1]]
            im = ax[i,j].imshow(I, extent=extent, aspect="auto", origin="lower", vmin=1-dT, vmax=1+dT)
            f.close()

    for i,Teff in enumerate(Teff_list):
        ax[i,-1].yaxis.set_label_position("right")
        ax[i,-1].set_ylabel(r"$T_{eff}$ = "+f"{Teff} K")
        ax[i,0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax[i,0].set_ylabel(r"$y$ [Mm]")
    for j, B0 in enumerate(B0_list):
        ax[-1,j].set_xlabel(fr"$B_0$ = {B0} G")
    
    fig.colorbar(im, ax=ax.ravel().tolist(), location="top", fraction=0.046, pad=0.02, label=label)
    plt.subplots_adjust(top=0.84, hspace=0, wspace=0)

    plt.savefig("figures/int_imgs/"+figname, bbox_inches="tight")
    
def plot_Tb_imgs(wave, figname, label, dT):

    fs = np.array([15,15])*0.5
    fig,ax = plt.subplots(5,5,figsize=(fs[0],fs[1]))

    Teff_list = [6500,5770,5000,4000,3200]
    B0_list = [0,100,200,500]

    for i,Teff in enumerate(Teff_list):
        ax[i,-1].get_xaxis().set_ticks([])
        ax[i,-1].yaxis.tick_right()
        ax[i,-1].yaxis.set_label_position("right")
        #ax[i,-1].set_xlim([0,1.1])
        color = f"C{i}"
        
        #F = h5py.File(input_files[Teff][0], "r")
        #nx = 1/(np.array(F["dx"][0])*1e-8)
        
        for j, B0 in enumerate(B0_list):
            f = h5py.File(output_files[Teff][B0][wave], "r")
            F = h5py.File(input_files[Teff][B0], "r")
            x = np.cumsum(np.array(F["dx"]))*1e-8

            intensity = ((wave*u.mm)**2/c)*np.array(f["intensity"])*I_unitsw
            
            Tb = brightness_temperature(intensity, wave*u.mm).to(u.kK)

            #I = intensity/np.mean(intensity)

            ax[i,j].get_xaxis().set_ticks([])
            ax[i,j].get_yaxis().set_ticks([])
            
            extent = [x[0],x[-1],x[0],x[-1]]
            im = ax[i,j].imshow(Tb.value/np.mean(Tb.value), extent=extent, aspect="auto", origin="lower", cmap="inferno",
                                vmin=1-dT, vmax=1+dT)
            
            if B0==0:
                x0 = 0.1*x[-1]
                y0 = 0.9*x[-1]
                ax[i,j].plot([x0, x0+1], [y0,y0], lw=5, color="black")
                ax[i,j].plot([x0, x0+1], [y0,y0], lw=3, color="white")
            
            if B0==0:
                Tmax = np.max(Tb.value)
                Tmin = np.min(Tb.value)
                DT = Tmax-Tmin
                Trange = (Tmin - DT/5, Tmax + DT/5)
                ax[i,-1].set_ylim([Trange[0], Trange[-1]])
            
            Tb_hist, bin_edges = np.histogram(Tb.value.flatten(), density=True, range=Trange, bins=30)
            
            #ax[i,-1].stairs(Tb_hist/np.max(Tb_hist), bin_edges, orientation="horisontal", color=color, ls=ls_list[j], lw=1)
            bin_centres = (bin_edges[1:] + bin_edges[:-1])/2
            #ax[i,-1].plot(Tb_hist/np.max(Tb_hist), bin_centres, color=color, ls=ls_list[j], lw=1)
            ax[i,-1].plot(Tb_hist, bin_centres, color=color, ls=ls_list[j], lw=1)
            f.close()
            
    
    
    for i,Teff in enumerate(Teff_list):
        #ax[i,-1].yaxis.set_label_position("right")
        ax[i,0].set_ylabel(r"$T_{eff}$ = "+f"{Teff} K")
        ax[i,-1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax[i,-1].set_xlim(left=0)
        #ax[i,0].set_ylabel(r"$y$ [Mm]")
    for j, B0 in enumerate(B0_list):
        ax[-1,j].set_xlabel(fr"$B_0$ = {B0} G")
        
    ax[-1,-1].set_xlabel(r"Histogram, $T_B$ [kK]")
    
    fig.colorbar(im, ax=ax[:,:-1].ravel().tolist(), location="top", fraction=0.046, pad=0.02, label=label)
    plt.subplots_adjust(top=0.84, hspace=0, wspace=0)

    plt.savefig("figures/Tb_imgs/"+figname, bbox_inches="tight")

#plot_int_imgs(0.0005, "int_500nm.pdf", r"$I_{500 nm}/\langle I_{500 nm} \rangle$", dT=0.35)
#plot_int_imgs(0.3, "int_00_30mm.pdf", r"$I_{0.3 mm}/\langle I_{0.3 mm} \rangle$", dT=0.15)
#plot_int_imgs(1, "int_01_00mm.pdf",   r"$I_{1.0 mm}/\langle I_{1.0 mm} \rangle$", dT=0.20)
#plot_int_imgs(3, "int_03_00mm.pdf",   r"$I_{3.0 mm}/\langle I_{3.0 mm} \rangle$", dT=0.25)
#plot_int_imgs(6, "int_06_00mm.pdf",   r"$I_{6.0 mm}/\langle I_{6.0 mm} \rangle$", dT=0.30)
#plot_int_imgs(9, "int_09_00mm.pdf",   r"$I_{9.0 mm}/\langle I_{9.0 mm} \rangle$", dT=0.35)

plot_Tb_imgs(0.0005, "Tb_500nm.pdf", r"$T_{B, 500 nm}/\langle T_{B, 500 nm} \rangle$",   dT=0.075)
plot_Tb_imgs(0.3,    "Tb_00_30mm.pdf", r"$T_{B, 0.3 mm}/\langle T_{B, 0.3 nm} \rangle$", dT=0.15)
plot_Tb_imgs(1,      "Tb_01_00mm.pdf", r"$T_{B, 1.0 mm}/\langle T_{B, 1.0 mm} \rangle$", dT=0.15)
plot_Tb_imgs(3,      "Tb_03_00mm.pdf", r"$T_{B, 3.0 mm}/\langle T_{B, 3.0 mm} \rangle$", dT=0.25)
plot_Tb_imgs(6,      "Tb_06_00mm.pdf", r"$T_{B, 6.0 mm}/\langle T_{B, 6.0 mm} \rangle$", dT=0.25)
plot_Tb_imgs(9,      "Tb_09_00mm.pdf", r"$T_{B, 9.0 mm}/\langle T_{B, 9.0 mm} \rangle$", dT=0.30)