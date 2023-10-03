import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import h5py
#import uio

from astropy import units as u
from astropy import constants

from scipy import interpolate
from scipy import stats

c = constants.c
kB = constants.k_B

art_output = h5py.File("/uio/hume/student-u19/jonasrth/ART1/example/solar_synthetic2.h5","r")
art_input = h5py.File("/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/d3t57g44c_v000G_n013_end_art.h5", "r")
#cobold_file = uio.File("/mn/stornext/d19/RoCS/svenwe/jonast/cobold/model/d3t57g44c_v000G_n013.end")

intensity = np.array(art_output["Stokes_I"][0,:,:,:])*u.erg*u.cm**(-2)*u.s**(-1)*u.Hz**(-1)
wave = np.array(art_output["Wavelength"])*u.AA
zt1 = np.array(art_output["Tau1"])[0,:,:,:]*u.cm

z = np.array(art_input["z"][0,0,0,:])*u.cm

shape = intensity.shape
dx = np.array(art_input["dx"][0])*1e-8
dy = dx

extent = [0,shape[0]*dx,0,shape[1]*dy]

def brightness_temp(intensity,wave):
    """
    intensity in cgs units (per frequency)
    wavelength in Ångström
    """
    Tb = (wave**2/(2*kB)*intensity).to("K")
    return Tb

fig,ax = plt.subplots(2,2,figsize=(7.5,7))
for i in range(4):
    x = brightness_temp(intensity[:,:,i],wave[i]).value.flatten()
    
    inds = np.argmin(abs(zt1[:,:,i][:,:,np.newaxis] - z[np.newaxis,np.newaxis,:]), axis=2)
    
    Y = np.take_along_axis(np.array(art_input["temperature"])[0,...], inds[...,np.newaxis], axis=2)[...,0]
    y = Y.flatten()
    
    res = stats.linregress(x, y)
    
    t = np.linspace(0,10000, 100)
    label=fr"$({np.round(res.slope,2)}\pm{np.round(res.stderr,3)})T_b + ({int(res.intercept)}\pm{int(res.intercept_stderr)})$ K"

    H, x_edges, y_edges = np.histogram2d(x,y,bins=[75,75])
    ax.flat[i].plot(t, res.slope*t + res.intercept, ls="dashed", color="red", label=label)
    #ax.flat[i].plot(x_edges, 1*x_edges + 0, ls="dashed", color="green")

    extent=[np.min(x),np.max(x),np.min(y),np.max(y)]
    #ax.flat[i].imshow(np.log10(H.T + 1), origin="lower", extent=extent, aspect="auto", cmap="Greys")
    #ax.flat[i].contourf([x_edges[:-1], y_edges[:-1]], np.log10(H.T+1), extent=extent, cmap = "Greys")
    ax.flat[i].contourf(np.log10(H.T+1), extent=extent, cmap = "Greys", extend="min")
    #ax.flat[i].imshow(H.T, origin="lower", extent=extent, aspect="auto", cmap="Greys")
    #ax.flat[i].set_xlabel(r"$T_b$ [K]")
    #ax.flat[i].set_ylabel(r"T [K]")
    #ax.flat[i].set_title(r"$\lambda =$ " + f"{np.round(wave[i].to(u.mm),1)}")
    ax.flat[i].set_ylim([np.min(extent), np.max(extent)])
    ax.flat[i].set_xlim([np.min(extent), np.max(extent)])
    ax.flat[i].text(0.55,0.1, r"$\lambda =$ " + f"{np.round(wave[i].to(u.mm),1)}", size=12,transform=ax.flat[i].transAxes)
    ax.flat[i].legend(loc="upper center")

for i in range(2):
    #ax[-1,i].set_xlabel(r"$T_b$ [K]")
    ax[-1,i].set_xlabel("Brightness temperature [K]")
    ax[i,0].set_ylabel("Gas temperature [K]") 
    
plt.tight_layout()

figname = "TbvsTg_contour.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")