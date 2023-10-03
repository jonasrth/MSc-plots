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


data_path = "/mn/stornext/d20/RoCS/atulm/Project1_stars/SED_data/"
SED_path = data_path + "Clean_SED_data/"

path_art = "/mn/stornext/d19/RoCS/jonasrth/ART/SED/"

def get_star_data(star_name):
    """
    Collects necessary data on the stars to compare them with model SED.
    star_name (str) must fit with one of the 12 stars compiled in the EMISSA project.
    """
    
    # Collecting SED data for star_name
    filename = star_name + "_CleanSED.csv"
    with open(SED_path+filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        csv_array = np.array(list(csv_reader))
        
        n_freq = np.argwhere(csv_array[0]=="sed_freq")[0][0]
        n_flux = np.argwhere(csv_array[0]=="sed_flux")[0][0]
        n_eflux = np.argwhere(csv_array[0]=="sed_eflux")[0][0]
        n_tab = np.argwhere(csv_array[0]=="_tabname")[0][0]

        sed_freq = csv_array[1:,n_freq].astype(np.float64)*u.GHz
        sed_flux = csv_array[1:,n_flux].astype(np.float64)*u.Jy
        sed_eflux = csv_array[1:,n_eflux].astype(np.float64)*u.Jy
        tabname = csv_array[1:,n_tab]
    
    # Collecting radius and distance data for star_name
    with open(data_path + "star_props.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        csv_array = np.array(list(csv_reader))

        n_d = np.argwhere(csv_array[0]==" Distance (pc)")[0][0]
        n_ed = np.argwhere(csv_array[0]==" Dist_err")[0][0]
        n_R = np.argwhere(csv_array[0]=="Radius (Rs)")[0][0]
        n_eR = np.argwhere(csv_array[0]=="Rad_err")[0][0]

        m = np.argwhere(csv_array[:,1]==star_name.replace("_"," "))[0][0]

        d = float(csv_array[m, n_d])*u.pc
        d_err = float(csv_array[m, n_ed])*u.pc
        R = float(csv_array[m, n_R])*const.R_sun
        R_err = float(csv_array[m, n_eR])*const.R_sun
    
    # Returning collected data in dictionary:
    
    data = {}
    
    data["sed_freq"]  = sed_freq
    data["sed_flux"]  = sed_flux
    data["sed_eflux"] = sed_eflux
    data["tabname"] = tabname
    
    data["d"]     = d
    data["d_err"] = d_err
    data["R"]     = R
    data["R_err"] = R_err
    
    return data  

def subplot_SED(star_name, model_data, gs):
    """
    Code for each subplot.
    star_name (str) must fit with one of the 12 stars compiled in the EMISSA project. 
    model_data (dict={int : .h5 file}) includes SED data for model atmospheres
    gs (class) specifies axes in subplot
    """
    
    data = get_star_data(star_name)
    
    mod_int = np.array([])
    mod_wav = np.array([])
    for file in model_data.values():
        mod_int = np.append(mod_int, np.mean(np.array(file["Stokes_I"][0,...]),axis=(0,1)))
        mod_wav = np.append(mod_wav, np.array(file["Wavelength"]))
    
    mod_freq = (c/(mod_wav*u.angstrom)).to(u.GHz)
    
    mod_freq, inds = np.unique(mod_freq, return_index=True)
    mod_int = mod_int[inds]*I_units
    
    mod_flux = (np.pi*(data["R"]/data["d"])**2*mod_int*u.sr).to(u.Jy)
    mod_flux_max = (np.pi*((data["R"]+data["R_err"])/(data["d"]-data["d_err"]))**2*mod_int*u.sr).to(u.Jy)
    mod_flux_min = (np.pi*((data["R"]-data["R_err"])/(data["d"]+data["d_err"]))**2*mod_int*u.sr).to(u.Jy)
    
    a = np.argmin(abs(mod_freq.value - np.min(data["sed_freq"]).value)) - 1
    b = np.argmin(abs(mod_freq.value - np.max(data["sed_freq"]).value)) + 1
    
    ### Interpolation
    
    f = interp1d(mod_freq, mod_flux)
    
    ### Plotting:
    
    #plt.suptitle(star_name.replace("_"," "), x=0.5, y=1.0)
    
    ax0 = fig.add_subplot(gs[1])
    
    ax0.text(0.1, 0.9, star_name.replace("_"," "), transform=ax0.transAxes)
    
    cmap = cm.get_cmap("gnuplot")
    
    it = 0
    gradient = np.linspace(0,1, len(unique(data["tabname"])))
    for tab in unique(data["tabname"]):
        
        n = np.argwhere(data["tabname"]==tab)[:,0]
        ax0.errorbar(data["sed_freq"][n].value, data["sed_flux"][n].value, yerr=data["sed_eflux"][n].value,
                     color=cmap(gradient[it]), ls="None", marker="o", label=tab)
        it += 1
    
    handles, labels = ax0.get_legend_handles_labels()
    legd=ax0.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5,1.05), ncol=4)
    
    ax0.plot(mod_freq[a:b], mod_flux[a:b], color="black", ls="solid")
    ax0.fill_between(mod_freq[a:b], y1=mod_flux_max[a:b], y2=mod_flux_min[a:b], color="grey", alpha=0.5)
    ax0.axvspan(0,1e3, color="grey", alpha=0.2)
    ax0.set_ylabel("Flux [Jy]")
    
    ax0.xaxis.grid(which="both")
    ax0.yaxis.grid(which="major")
    ax0.set_yscale("log")
    
    ax1 = fig.add_subplot(gs[2], sharex=ax0)
    ax1.errorbar(data["sed_freq"], (data["sed_flux"].value - f(data["sed_freq"]))/f(data["sed_freq"]),
                color="black", ls="None", marker="o")
    ax1.axvspan(0,1e3, color="grey", alpha=0.2)
    ax1.set_ylabel(r"$\Delta S/S_{mod}$")
    ax1.set_xlabel("Frequency [GHz]")
    ax1.xaxis.grid(which="both")
    ax1.yaxis.grid(which="major")
    
    ax1.set_xscale("log")
    plt.setp(ax0.get_xticklabels(),visible=False)
    #plt.subplots_adjust(hspace=0)
    
star_name_list = ["Gam_Lep", "Alf_Cen_A", "Alf_Cen_B", "GJ_2006_A"]
model_data_list = [{0 : h5py.File(path_art + "d3t65g45_000G_SED.h5","r")},
                   {0 : h5py.File(path_art + "d3t57g44_000G_SED.h5","r")},
                   {0 : h5py.File(path_art + "d3t50g45_000G_SED.h5","r")},
                   {0 : h5py.File(path_art + "d3t32g45_000G_SED.h5","r")}]

fig = plt.figure(figsize=(20,20))

outer = gridspec.GridSpec(2, 2, wspace=0.2, hspace=0.2) 

for i in range(4):
    inner = gridspec.GridSpecFromSubplotSpec(3, 1, hspace=0, height_ratios=[1, 4, 1],
                                            subplot_spec=outer[i])
    
    subplot_SED(star_name_list[i], model_data_list[i], gs=inner)
    

figname = "SED_EMISSA_plot.pdf"   
plt.savefig("figures/" + figname, bbox_inches="tight")