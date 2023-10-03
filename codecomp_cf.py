import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import cumtrapz

from astropy import units as u
from astropy import constants as const

c = const.c

import h5py

import sys

sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

TR_dict = np.load("files/TR_dict.pkl", allow_pickle=True)

#Lf3D_wave_list = ["5000", "00_30", "01_00", "03_00", "06_00", "09_00"]
art_wave_list = ["500nm", "00_30mm", "01_00mm", "03_00mm", "06_00mm", "09_00mm"]
wave_list = [0.0005, 0.3, 1, 3, 6, 9]

art_path = "/mn/stornext/d19/RoCS/jonasrth/ART/code_comparison/"
#rh_path = "/mn/stornext/d19/RoCS/jonasrth/RH_input/output/d3t57g44c_v000G_n013_it000_full_rh/output/"
rh_path = "/mn/stornext/d19/RoCS/jonasrth/RH/code_comparison/d3t57g44c_v000G_n013_it000/"
#rh_path = "/uio/hume/student-u19/jonasrth/rh/rh15d/run_compare/output/"


F = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/d3t57g44c_v000G_n013_it000_full_art.h5", "r")
z = np.array(F["z"][0,0,0,:])*u.cm
logTR = TR_dict[5700][0]
ks = np.argmin(abs(logTR))
kc = np.argmin(abs(logTR+5))
    
z = z - z[ks]

I_units_rh = u.W*u.m**(-2)*u.Hz**(-1)*u.sr**(-1)
I_units_art = u.erg*u.cm**(-2)*u.s**(-1)*u.Hz**(-1)*u.sr**(-1)
I_units_Lf = u.erg*u.cm**(-2)*u.s**(-1)*u.angstrom**(-1)*u.sr**(-1)



def compare_CFs(wave_index, figname):
    """
    """
    
    i = wave_index
    
    # Output files:
    f_L = h5py.File(output_files[5700][0][wave_list[i]], "r")
    f_A = h5py.File(art_path + "d3t57g44_000G_CF_" + art_wave_list[i] + ".h5", "r")
    f_R = h5py.File(rh_path + "output_ray.hdf5", "r")
    
    d_CF = {}
    
    ### ART CF:
    CF = np.mean(np.array(f_A["CFunc"][0,...]), axis=(0,1))
    #CF = np.array(f_A["CFunc"][0,0,0,:])
    d_CF["A"] = CF*I_units_art/u.cm
    
    ### RH CF:
    chi = np.mean(np.array(f_R["chi"][...,i]), axis=(0,1))*u.m**(-1)
    Snu = np.mean(np.array(f_R["source_function"][...,i]), axis=(0,1))*I_units_rh

    tau = cumtrapz(chi.value, z.to(u.m).value, initial=0)
    integrand = chi*Snu*np.exp(tau)
    max_ = np.max(integrand, axis=0)
    d_CF["R"] = (integrand).to(I_units_art/u.cm)
    
    ### Linfor3D CF:
    z_CF = np.array(f_L["contfctz"])*u.cm
    a = np.argmax(z_CF)
    z_CF = z_CF[:a]
    CF = np.mean(np.array(f_L["contfunc"][:a]), axis=(1,2))
    #CF = np.array(f_L["contfunc"][:a,0,0])[::-1]
    d_CF["L"] = CF*I_units_Lf/u.cm
    d_CF["L"] = (((wave_list[i]*u.mm)**2/c)*d_CF["L"]).to(I_units_art/u.cm)
    
    fs = np.array([10,10])*0.5
    fig,ax = plt.subplots(figsize=(fs[0],fs[1]))
    
    ax.plot(z_CF.to(u.Mm), d_CF["L"]/np.max(d_CF["L"]), label="Linfor3D")
    ax.plot(z.to(u.Mm), d_CF["A"]/np.max(d_CF["A"]), label="ART")
    ax.plot(z.to(u.Mm), d_CF["R"]/np.max(d_CF["R"]), label="RH")
    #ax.title("Wavelength: " + str(wave_list[i]) + " mm")
    ax.set_xlabel("z [Mm]")
    ax.set_ylabel("Normalised contribution function")
    ax.set_xlim([0,2])
    ax.axvline(x=z[kc].to(u.Mm).value, color="black", ls="dashdot", lw=1.2)
    ax.legend()
    
    w = np.array(f_A["Wavelength"])*u.angstrom
    
    f_A.close()
    f_R.close()
    f_L.close()
    
    plt.savefig("figures/"+figname, bbox_inches="tight")
    
compare_CFs(3, figname="codecomp_cf.pdf") 