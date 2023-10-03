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

def plot_SED(star_name, model_data, figname):
    """
    """
    
    data = get_star_data(star_name)
    
    f_inds = np.argwhere(data["sed_freq"].value > 1e6)[:,0]
    data["sed_freq"] = np.delete(data["sed_freq"], f_inds)
    data["sed_flux"] = np.delete(data["sed_flux"], f_inds)
    
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
    #ax0.text(0.1, 0.9, star_name.replace("_"," "), transform=ax0.transAxes)
    
    ## Observed data:
    cmap = cm.get_cmap("gnuplot")
    gradient = np.linspace(0, 1, len(unique(data["tabname"])))
    data["tabname"] = np.delete(data["tabname"], f_inds)
    for it,tab in enumerate(unique(data["tabname"])):
        
        n = np.argwhere(data["tabname"]==tab)[:,0]
        ax0.errorbar(data["sed_freq"][n].value, data["sed_flux"][n].value, yerr=data["sed_eflux"][n].value,
                     color=cmap(gradient[it]), ls="None", marker="o", label=tab)
        if tab=="ALMA data":
            n_ALMA = it
    
    ## Model data:
    for B0 in model_data.keys():
        
        i_min = np.argmin(abs(mod_freq[B0].value - np.min(data["sed_freq"]).value)) - 1
        i_max = np.argmin(abs(mod_freq[B0].value - np.max(data["sed_freq"]).value)) + 2
        
        i_min = np.argmin(abs(mod_freq[B0].value - 7)) - 1
        
        label = fr"$B_0$ = {B0} G"
        ax0.plot(mod_freq[B0][i_min:i_max], mod_flux[B0][i_min:i_max], ls="solid", label=label)
        #ax0.fill_between(mod_freq[a:b], y1=mod_flux_max[a:b], y2=mod_flux_min[a:b], color="grey", alpha=0.5)
    
    handles, labels = ax0.get_legend_handles_labels()
    data_len = len(model_data.keys())
    handle_list = [handles[n_ALMA+data_len]]
    label_list = [labels[n_ALMA+data_len]]
    for i in range(data_len):
        handle_list.append(handles[i])
        label_list.append(labels[i])
        
    #handle_list = [handles[n_ALMA+4], handles[0], handles[1], handles[2], handles[3]]
    #label_list = [labels[n_ALMA+4], labels[0], labels[1], labels[2], labels[3]]
    legd=ax0.legend(handle_list, label_list, loc="lower center", bbox_to_anchor=(0.5,1.01), ncol=5)
    
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
        ax1.errorbar(data["sed_freq"], rel_error, ls="None", marker="o")
    
    ax1.axvspan(0,1e3, color="grey", alpha=0.2)
    ax1.set_ylabel(r"$\Delta S/S_{mod}$")
    ax1.set_xlabel("Frequency [GHz]")
    ax1.xaxis.grid(which="both")
    ax1.yaxis.grid(which="major")
    
    ax1.set_xscale("log")
    plt.setp(ax0.get_xticklabels(),visible=False)
    
    plt.savefig("figures/" + figname, bbox_inches="tight")
    
    
    

star_name_list = ["Gam_Vir_A", "Gam_Vir_B", "Eta_Crv", "Gam_Lep", "Alf_Cen_A", "61_Vir", "Alf_Cen_B", "Eps_Eri", "GJ_2006_A", "Proxima_Cen"]
star_letter_list = ["C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]

model_name_list = ["t65", "t65", "t65", "t65", "t57", "t57", "t50", "t50", "t32", "t32"]

for i in range(len(star_name_list)):
    print(star_letter_list[i] + ": " +star_name_list[i]+" - "+model_name_list[i])

SEDs = {"t65" : {  0 : {0 : path_art + "d3t65g45_000G_SED.h5"},
                 100 : {0 : path_art + "d3t65g45_100G_SED.h5"},
                 200 : {0 : path_art + "d3t65g45_200G_SED.h5"},
                 500 : {0 : path_art + "d3t65g45_500G_SED.h5"},
                },
        "t57" : {  0 : {0 : path_art + "d3t57g44_000G_SED.h5"},
                 100 : {0 : path_art + "d3t57g44_100G_SED.h5"},
                 200 : {0 : path_art + "d3t57g44_200G_SED.h5"},
                 500 : {0 : path_art + "d3t57g44_500G_SED.h5"},
                },
        "t50" : {  0 : {0 : path_art + "d3t50g45_000G_SED.h5"},
                 100 : {0 : path_art + "d3t50g45_100G_SED.h5"},
                 200 : {0 : path_art + "d3t50g45_200G_SED.h5"},
                 500 : {0 : path_art + "d3t50g45_500G_SED.h5"},
                },
        "t32" : {  0 : {0 : path_art + "d3t32g45_000G_SED.h5"},
                 100 : {0 : path_art + "d3t32g45_100G_SED.h5"},
                 200 : {0 : path_art + "d3t32g45_200G_SED.h5"},
                 500 : {0 : path_art + "d3t32g45_500G_SED.h5"},
                }
       }
"""
for i in range(len(star_name_list)):
    figname = "EMISSA/SED_"+star_letter_list[i]+"_"+model_name_list[i]+".pdf"
    plot_SED(star_name_list[i], SEDs[model_name_list[i]], figname)
"""   
for i in range(len(star_name_list)):
    figname = "presentation/SED_"+star_letter_list[i]+"_"+model_name_list[i]+".pdf"
    plot_SED(star_name_list[i], SEDs[model_name_list[i]], figname)