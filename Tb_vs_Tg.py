import matplotlib.pyplot as plt
import numpy as np

import h5py

from astropy import units as u
from astropy import constants

I_units = u.erg*u.cm**(-2)*u.s**(-1)*u.Hz**(-1)*u.sr**(-1)

from scipy import interpolate
from scipy import stats

c = constants.c
kB = constants.k_B
h = constants.h

### Output and input files:

"""
f = {
    "500nm" : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/test_mode1_500nm.h5", "r"),
    "1mm"   : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/test_mode1_1mm.h5", "r"),
    "3mm"   : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/test_mode1_3mm.h5", "r")
    }
"""

f = {
    "500nm" : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/new_mode1_500nm.h5", "r"),
    "1mm"   : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/new_mode1_1mm.h5", "r"),
    "3mm"   : h5py.File("/mn/stornext/d19/RoCS/jonasrth/ART/test_modes_t57/new_mode1_3mm.h5", "r")
    }

path_Lf3D = "/mn/stornext/d19/RoCS/alma/emissa_sim/linfor3D/outhdf/"
f_Lf3D = {"500nm" : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_05000_mu1_00_linfor_3D_2.hdf", "r"),
          "1mm"   : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_01mm_mu1_00_linfor_3D_2.hdf", "r"),
          "3mm"   : h5py.File(path_Lf3D + "d3t57g44c_v000G_n019_it000_03mm_mu1_00_linfor_3D_2.hdf", "r")}

#F = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/d3t57g44c_v000G_n019_it000_full_art.h5", "r")
F = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/tst/d3t57g44_v000G_n019_art_it000_mode1.h5", "r")


### Extracting quantities:

intensity, wave, zt1 = {}, {}, {}
for key in f.keys():
    intensity[key] = np.array(f[key]["Stokes_I"][0,:,:,:])*I_units
    wave[key] = np.array(f[key]["Wavelength"][0])*u.AA
    zt1[key] = np.array(f[key]["Tau1"])[0,:,:,0]*u.cm
    
    """
    z_CF = np.array(f_Lf3D[key]["contfctz"])*1e-8
    k = np.argmax(z_CF)
    z_CF = z_CF[:k+1]

    CF = np.array(f_Lf3D[key]["contfunc"][:k+1,:,:])
    
    zt1[key] = np.trapz(z_CF[:,np.newaxis,np.newaxis]*CF, z_CF, axis=0)/np.trapz(CF, z_CF, axis=0)
    """

z = np.array(F["z"][0,0,0,:])*u.cm
shape = intensity[key].shape
dx = np.array(F["dx"][0])*1e-8
dy = dx

extent = [0,shape[0]*dx,0,shape[1]*dy]

def brightness_temperature(intensity,wave):
    """
    intensity in cgs units (per frequency)
    wavelength in Ångström
    """
    #Tb = (wave**2/(2*kB)*intensity*u.sr).to("K")
    Tb = (h*c/(wave*kB)*np.log(1 + (2*h*c/(wave**3*intensity*u.sr)).to(1))**(-1)).to("K")
    return Tb

fig,ax = plt.subplots(1,3,figsize=(12,4))
for i,key in enumerate(f.keys()):
    
    x = brightness_temperature(intensity[key],wave[key]).value.flatten()
    
    inds = np.argmin(abs(zt1[key][:,:,np.newaxis] - z[np.newaxis,np.newaxis,:]), axis=2)
    
    Y = np.take_along_axis(np.array(F["temperature"])[0,...], inds[...,np.newaxis], axis=2)[...,0]
    y = Y.flatten()
    
    res = stats.linregress(x, y)
    
    t = np.linspace(0,10000, 100)
    ax[i].plot(t, res.slope*t + res.intercept, ls="dashed", color="grey", 
               label=fr"$({np.round(res.slope,2)}\pm{np.round(res.stderr,3)})T_b + ({int(res.intercept)}\pm{int(res.intercept_stderr)})$ K")
    
    #print(np.sign(res.intercept))

    #H = np.histogram2d(x,y,bins=[50,50])[0]

    extent=[np.min(x),np.max(x),np.min(y),np.max(y)]
    #ax.flat[i].imshow(np.log10(H.T + 1), origin="lower", extent=extent, aspect="auto")
    ax[i].scatter(x, y, color="black", alpha=1/510)
    #ax.flat[i].imshow(H.T, origin="lower", extent=extent, aspect="auto")
    
    dx = extent[1] - extent[0]
    dy = extent[3] - extent[2]
    ax[i].set_xlim(extent[0] - dx/20, extent[1] + dx/20)
    ax[i].set_ylim(extent[2] - dy/20, extent[3] + dy/20)
    
    ax[i].set_axisbelow(True)
    ax[i].grid()
    ax[i].set_xlabel(r"Brightness temperature [K]")
    ax[i].set_title(r"$\lambda =$ " + f"{key[:-2]} {key[-2:]}")
    
    ax[i].legend(loc="upper left")

ax[0].set_ylabel(r"Gas temperature at $z(\tau=1)$ [K]")
plt.tight_layout()

figname = "Tb_vs_Tg.png"
plt.savefig("figures/" + figname, bbox_inches="tight")