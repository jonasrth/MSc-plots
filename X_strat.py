import numpy as np
import matplotlib.pyplot as plt

import pickle
import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

HP_dict = np.load("files/HP_dict.pkl", allow_pickle=True)
TR_dict = np.load("files/TR_dict.pkl", allow_pickle=True)

def plot_T_strat(figname):
    """
    """
    
    fig,ax = plt.subplots(figsize=(10,8))

    ls_list = ["solid", "dashed", "dashdot", "dotted"]
    for i, Teff in enumerate([6500,5770,5000,4000,3200]):
        color = f"C{i}"
        for j,B0 in enumerate([0,100,200,500]):
            F = h5py.File(input_files[Teff][B0], "r")
            
            T = np.mean(np.array(F["temperature"][0,...]), axis=(0,1))
            #logP = HP_dict[Teff][B0]
            logTR = TR_dict[Teff][B0]

            if B0==0:
                # For labelling
                ax.plot(logTR, T*1e-3, color=color, ls=ls_list[j], label=r"$T_{eff}$ = "f"{Teff} K")
                #ax.plot(logP, T*1e-3, color=color, ls=ls_list[j], label=f"{Teff} K")
            ax.plot(logTR, T*1e-3, color=color, ls=ls_list[j])
            #ax.plot(logP, T*1e-3, color=color, ls=ls_list[j])

    ax.legend()
    ax.axvline(x = 0.0, color="k", ls="dashed", lw=1.2)
    ax.axvline(x = -5.0, color="k", ls="dashdot", lw=1.2)
    ax.set_xlabel(r"$\log_{10}(\langle \tau_R \rangle_z)$")
    #ax.set_xlabel(r"$\langle P \rangle_z/P_0$")
    ax.set_ylabel("Gas temperature [kK]")
    ax.grid()
    
    plt.savefig("figures/"+figname, bbox_inches="tight")
    
def plot_Tstd_strat(figname):
    """
    """
    
    fig,ax = plt.subplots(figsize=(10,8))

    ls_list = ["solid", "dashed", "dashdot", "dotted"]
    for i, Teff in enumerate([6500,5770,5000,4000,3200]):
        color = f"C{i}"
        for j,B0 in enumerate([0,100,200,500]):
            F = h5py.File(input_files[Teff][B0], "r")
            
            T = np.mean(np.array(F["temperature"][0,...]), axis=(0,1))
            
            T_rms = np.sqrt(np.mean(((np.array(F["temperature"][0,...]) - T)/T)**2, axis=(0,1)))
            
            T_std = np.std(np.array(F["temperature"][0,...]), axis=(0,1))
            T_ = T_std/T
            #T_ = T_rms
            
            #logP = HP_dict[Teff][B0]
            logTR = TR_dict[Teff][B0]

            if B0==0:
                # For labelling
                ax.plot(logTR, T_, color=color, ls=ls_list[j], label=r"$T_{eff}$ = "f"{Teff} K")
                #ax.plot(logP, T_, color=color, ls=ls_list[j], label=f"{Teff} K")
            ax.plot(logTR, T_, color=color, ls=ls_list[j])
            #ax.plot(logP, T_, color=color, ls=ls_list[j])  

    ax.legend()
    ax.axvline(x = 0, color="k", ls="dashed", lw=1.2)
    ax.axvline(x = -5.0, color="k", ls="dashdot", lw=1.2)
    ax.set_ylabel(r"$\sqrt{\langle (T - \langle T \rangle_z)^2\rangle_z}/\langle T \rangle_z$")
    ax.set_xlabel(r"$\log_{10}(\langle \tau_R \rangle_z)$")
    #ax.set_xlabel(r"$\langle P \rangle_z/P_0$")
    ax.set_ylim([0,0.6])
    ax.grid()
    
    plt.savefig("figures/"+figname, bbox_inches="tight")
    
def plot_B_strat(figname):
    """
    """
    
    fig,ax = plt.subplots(figsize=(10,8))

    ls_list = ["solid", "dashed", "dashdot", "dotted"]
    for i, Teff in enumerate([6500,5770,5000,4000,3200]):
        color = f"C{i}"
        for j,B0 in enumerate([0,100,200,500]):
            F = h5py.File(input_files[Teff][B0], "r")
            
            B = np.mean(np.sqrt(np.array(F["bx"][0,...])**2 + np.array(F["by"][0,...])**2 + np.array(F["bz"][0,...])**2), axis=(0,1))
            
            #logP = HP_dict[Teff][B0]
            logTR = TR_dict[Teff][B0]

            if B0==0:
                # For labelling
                #ax.plot(logTR, B-B0, color=color, ls=ls_list[j], label=f"{Teff} K")
                ax.plot(logTR, B, color=color, ls=ls_list[j], label=r"$T_{eff}$ = "f"{Teff} K")
            #ax.plot(logTR, B-B0, color=color, ls=ls_list[j]) 
            ax.plot(logTR, B, color=color, ls=ls_list[j])

    ax.legend()
    ax.axvline(x = 0, color="k", ls="dashed", lw=1.2)
    ax.axvline(x = -5.0, color="k", ls="dashdot", lw=1.2)
    ax.set_ylabel("Magnetic field strength [G]")
    ax.set_xlabel(r"$\log_{10}(\langle \tau_R \rangle_z)$")
    #ax.set_xlabel(r"$\langle P \rangle_z/P_0$")
    ax.grid()
    
    plt.savefig("figures/"+figname, bbox_inches="tight")
    
def plot_z_strat(figname):
    """
    """
    
    fig,ax = plt.subplots(figsize=(10,8))

    ls_list = ["solid", "dashed", "dashdot", "dotted"]
    for i, Teff in enumerate([6500,5700,5000,4000,3200]):
        color = f"C{i}"
        for j,B0 in enumerate([0,100,200,500]):
            F = h5py.File(input_files[Teff][B0], "r")
            
            z = np.array(F["z"][0,0,0,:])*1e-8
            
            logTR = TR_dict[Teff][B0]

            if B0==0:
                # For labelling
                ax.plot(logTR, z, color=color, ls=ls_list[j], label=f"{Teff} K")
                #ax.plot(logP, B, color=color, ls=ls_list[j], label=f"{Teff} K")
            ax.plot(logTR, z, color=color, ls=ls_list[j]) 
            #ax.plot(logP, B, color=color, ls=ls_list[j)

    ax.legend()
    ax.axvline(x = 0, color="k", ls="dashed", lw=1.2)
    ax.set_ylabel("Vertical extent [Mm]")
    ax.set_xlabel(r"$\log_{10}(\langle \tau_R \rangle_z)$")
    #ax.set_xlabel(r"$\langle P \rangle_z/P_0$")
    ax.grid()
    
    plt.savefig("figures/"+figname, bbox_inches="tight")
        
    
plot_T_strat("T_strat.pdf")
plot_Tstd_strat("Tstd_strat.pdf")
plot_B_strat("B_strat.pdf")
#plot_z_strat("z_strat.pdf")