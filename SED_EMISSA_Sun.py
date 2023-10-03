import numpy as np

from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import gridspec

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from pandas import unique

import csv
import h5py

from astropy import constants as const
from astropy import units as u

I_units = u.erg*u.cm**(-2)*u.s**(-1)*u.Hz**(-1)*u.sr**(-1)

h = const.h
c = const.c
kB = const.k_B

import sys

sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

data_path = "/mn/stornext/d20/RoCS/atulm/Project1_stars/SED_data/Thekaekara_SunasStar_spec.tsv"

path_art = "/mn/stornext/d19/RoCS/jonasrth/ART/SED/"


def plot_SED(model_data, figname):
    """
    """
    fact = 1e26
    
    ws, Fs, dFs =  np.genfromtxt(data_path, dtype=float, usecols=(0,1,2), unpack=True, delimiter="\t", skip_header=1)
    
    frqs = c.value/ws
    
    # Scaling the SED data which is in Wm^-2 nm^-1 unit using fact that dlamda=dnu nu^2/c
    Fs[:-6]=Fs[:-6]*c/frqs[:-6]**2*10**-9 
    dFs[:-6]=dFs[:-6]*c/frqs[:-6]**2*10**-9
    Fs *= fact
    dFs *= fact
    
    data = {}
    data["sed_freq"] = frqs*u.GHz
    data["sed_flux"] = Fs*u.Jy
    data["sed_eflux"] = dFs*u.Jy
    data["R"] = const.R_sun
    data["d"] = const.au
    
    f_inds = np.argwhere(data["sed_freq"].value > 1e6)[:,0]
    data["sed_freq"] = np.delete(data["sed_freq"], f_inds)
    data["sed_flux"] = np.delete(data["sed_flux"], f_inds)
    data["sed_eflux"] = np.delete(data["sed_eflux"], f_inds)
    
    mod_int = {}
    mod_wav = {}
    
    mod_freq = {}
    mod_flux = {}
    
    for B0 in model_data.keys():
        mod_int[B0] = np.array([])
        mod_wav[B0] = np.array([])
        for output_file in model_data[B0].values():
            f = h5py.File(output_file, "r")
            mod_int[B0] = np.append(mod_int[B0], np.mean(np.array(f["Stokes_I"][0,...]),axis=(0,1)))
            mod_wav[B0] = np.append(mod_wav[B0], np.array(f["Wavelength"]))
        
        mod_freq[B0] = (c/(mod_wav[B0]*u.angstrom)).to(u.GHz)
        
        # In case the same frequency appears twice in output data:
        mod_freq[B0], inds = np.unique(mod_freq[B0], return_index=True)
        mod_int[B0] = mod_int[B0][inds]*I_units
        
        mod_flux[B0] = (np.pi*(data["R"]/data["d"])**2*mod_int[B0]*u.sr).to(u.Jy)
        
    
    ### Plotting:
    
    fig = plt.figure(figsize=(8,6.4))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)
    
    ax0 = fig.add_subplot(gs[0])
    #ax0.text(0.1, 0.9, "Sun", transform=ax0.transAxes)
            
    ax0.errorbar(data["sed_freq"][:-6].value, data["sed_flux"][:-6].value, yerr=data["sed_eflux"][:-6].value, 
                 color="b", ls="None", marker="o", label="(1)")
    ax0.errorbar(data["sed_freq"][-6:-4].value, data["sed_flux"][-6:-4].value, yerr=data["sed_eflux"][-6:-4].value, 
                 color="k", ls="None", marker="o", label="(2)")
    ax0.errorbar(data["sed_freq"][-4:].value, data["sed_flux"][-4:].value, yerr=data["sed_eflux"][-4:].value, 
                 color="g", ls="None", marker="o", label="(3)")
    
    ## Model data:
    for B0 in model_data.keys():
        
        i_min = np.argmin(abs(mod_freq[B0].value - np.min(data["sed_freq"]).value)) - 1
        i_max = np.argmin(abs(mod_freq[B0].value - np.max(data["sed_freq"]).value)) + 2
        
        i_min = np.argmin(abs(mod_freq[B0].value - 7)) - 1
        
        label = fr"$B_0$ = {B0} G"
        ax0.plot(mod_freq[B0][i_min:i_max], mod_flux[B0][i_min:i_max], ls="solid", label=label)
        #ax0.plot(mod_freq[B0][:i_max], mod_flux[B0][:i_max], ls="solid", label=label)
        #ax0.fill_between(mod_freq[a:b], y1=mod_flux_max[a:b], y2=mod_flux_min[a:b], color="grey", alpha=0.5)
    
        
    legd=ax0.legend(loc="lower center", bbox_to_anchor=(0.5,1.01), ncol=5)
    
    ax0.axvspan(0,1e3, color="grey", alpha=0.2)
    ax0.set_ylabel("Flux [Jy]")
    
    ax0.xaxis.grid(which="both")
    ax0.yaxis.grid(which="major")
    ax0.set_yscale("log")
    
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    
    #for B0 in model_data.keys():
    #    mod_intp = interp1d(mod_freq[B0], mod_flux[B0])
    #    rel_error = (data["sed_flux"].value - mod_intp(data["sed_freq"]))/mod_intp(data["sed_freq"])
    #    ax1.errorbar(data["sed_freq"], rel_error, ls="None", marker="o")
        
    for B0 in model_data.keys():
        mod_intp = interp1d(np.log10(mod_freq[B0].value), np.log10(mod_flux[B0].value))
        S_mod = 10**mod_intp(np.log10(data["sed_freq"].value))
        rel_error = (data["sed_flux"].value - S_mod)/S_mod
        ax1.errorbar(data["sed_freq"].value, rel_error, ls="None", marker="o")
    
    ax1.axvspan(0,1e3, color="grey", alpha=0.2)
    ax1.set_ylabel(r"$\Delta S/S_{mod}$")
    ax1.set_xlabel("Frequency [GHz]")
    ax1.xaxis.grid(which="both")
    ax1.yaxis.grid(which="major")
    
    ax1.set_xscale("log")
    plt.setp(ax0.get_xticklabels(),visible=False)
    
    #plt.savefig("figures/EMISSA/" + figname, bbox_inches="tight")
    plt.savefig("figures/presentation/" + figname, bbox_inches="tight")
    
model_data = {  0 : {0 : path_art + "d3t57g44_000G_SED.h5"},
              100 : {0 : path_art + "d3t57g44_100G_SED.h5"},
              200 : {0 : path_art + "d3t57g44_200G_SED.h5"},
              500 : {0 : path_art + "d3t57g44_500G_SED.h5"}
             }
#plot_SED(model_data, "SED_S_t57.pdf")
plot_SED(model_data, "SED_S_t57_test.pdf")
    
