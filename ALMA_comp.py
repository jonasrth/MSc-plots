import numpy as np
import matplotlib.pyplot as plt

import radio_beam as rb
from scipy.ndimage import convolve

from astropy import units as u
from astropy import constants as const

I_unitsw = u.erg*u.cm**(-2)*u.s**(-1)*u.angstrom**(-1)*u.sr**(-1)

c = const.c
h = const.h
kB = const.k_B
AU = const.au

import h5py
from astropy.io import fits

import sys

sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

def beam_kernel_calulator(bmaj_obs,bmin_obs,bpan_obs,ART_pxsz,size=199):
    """
    Calculate the beam array using the observed beam to be used for convolving the ART data
    """
    beam = rb.Beam(bmaj_obs*u.deg,bmin_obs*u.deg,bpan_obs*u.deg)
    #beam_kernel = np.asarray(beam.as_kernel(pixscale=ART_pxsz*u.deg, x_size=199, y_size=199))
    beam_kernel = np.asarray(beam.as_kernel(pixscale=ART_pxsz*u.deg, x_size=size, y_size=size))
    return beam_kernel

def brightness_temperature(intensity,wave):
    """
    intensity in cgs units (per frequency)
    wavelength in Ångström
    """
    Tb = (h*c/(wave*kB)*np.log(1 + (2*h*c/(wave**3*intensity*u.sr)).to(1))**(-1)).to("K")
    return Tb

fits_path = "/mn/stornext/d19/RoCS/svenwe/jonast/almaobs/"
filename = "solaralma.b3.fba.20161222_141931-150707.2016.1.00423.S.level4.k.fits"

wave = 3

Trange = (3500,9500)
Tbin_width = 100
Tbins = (Trange[-1] - Trange[0])//Tbin_width

vmin, vmax = 0.70,1.30

ff = fits.open(fits_path + filename)
hd = ff[0].header
ff.close()

fs = np.array([10,7.5])*0.8
fig = plt.figure(figsize=(fs[0],fs[1]))
gs = fig.add_gridspec(3,4)#, hspace=0, wspace=0)

ax0, ax1, ax2, ax3 = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1]), fig.add_subplot(gs[1,0]), fig.add_subplot(gs[1,1])
axs = [ax0,ax1,ax2,ax3]
ax = fig.add_subplot(gs[:2, 2:])

ax_ = fig.add_subplot(gs[-1,:])

letter_list = ["a","b","c","d","e","f"]

for k,B0 in enumerate([0,100,200,500]):
    
    axs[k].text(0.8,0.8, letter_list[k], color="white", size=20, transform=axs[k].transAxes)
    
    if k==0 or k==2:
        axs[k].set_ylabel("Y [arcsec]")
    else:
        axs[k].get_yaxis().set_ticks([])
        
    if k < 2:
        axs[k].set_xlabel("X [arcsec]")
        axs[k].xaxis.set_label_position("top")
        axs[k].xaxis.tick_top()
    else:
        axs[k].get_xaxis().set_ticks([])
        
    if k==2 or k==3:
        axs[k].get_xaxis().set_ticks([])
        
    #ax[i,-1].get_xaxis().set_ticks([])
    #ax[i,-1].yaxis.tick_right()
    #ax[i,-1].yaxis.set_label_position("right")
    
    F = h5py.File(input_files[5770][B0], "r")
    dx = np.array(F["dx"])[0]*u.cm
    DX = np.sum(np.array(F["dx"]))*u.cm
    ART_pxsz = (dx/AU).to(1)*180/np.pi
    
    f = h5py.File(output_files[5770][B0][3], "r")
    intensity = ((wave*u.mm)**2/c)*np.array(f["intensity"])*I_unitsw 
    Tb = brightness_temperature(intensity, wave*u.mm)

    # create the beam
    bmaj_obs, bmin_obs, bpan_obs = hd["BMAJ"], hd["BMIN"], hd["BPA"]
    
    #beam_kernel = beam_kernel_calulator(bmaj_obs,bmin_obs,bpan_obs,ART_pxsz, size=199)
    beam_kernel = beam_kernel_calulator(bmaj_obs,bmin_obs,bpan_obs,ART_pxsz, size=199)

    conv_results = convolve(Tb, beam_kernel, mode="wrap")
    conv_results_limres = conv_results[::16,::16]
    
    DEG = (DX/AU).to(1)*180/np.pi*3600
    extent = [-DEG/2, DEG/2, -DEG/2, DEG/2]
    im = axs[k].imshow(conv_results_limres/np.mean(conv_results_limres), origin="lower", cmap="inferno",
                       extent=extent, vmin=vmin, vmax=vmax)
    
    H,x_edges = np.histogram(conv_results.flatten(), bins=Tbins, range=Trange, density=True)
    #ax_.stairs(H, x_edges, label=f"{letter_list[k]}: $B_0$ = {B0} G")
    ax_.stairs(H*100, x_edges, label=f"{letter_list[k]}")
    #ax_.stairs(H, x_edges, label=f"{B0} G")
    
F.close()
f.close()

ff = fits.open(fits_path + filename)
hd = ff[0].header
img = np.array(ff[0].data)[0,0,0,:,:].T
shape = img.shape
ff.close()

#img_ = img.flatten()
#img_ = img_[~np.isnan(img_)]
#H,x_edges = np.histogram(img_, bins=Tbins, range=Trange, density=True)
#ax_.stairs(H, x_edges, label="e")
#ax_.set_xlabel("Brightness temperature [K]")
#ax_.legend()

dy = abs(hd["CDELT1"])*3600
dx = abs(hd["CDELT2"])*3600

#extent = [-dx/2*shape[1], dx/2*shape[1], -dy/2*shape[0], dy/2*shape[0]]
#im = ax.imshow(img/np.nanmean(img), origin="lower", cmap="inferno", extent=extent, vmin=vmin, vmax=vmax)
#ax.set_ylim([-DEG, DEG])
#ax.set_xlim([-DEG, DEG])

# resizing ALMA array to extent -DEG to +DEG in x and y

Npx = DEG/dx
Npy = DEG/dy

nx0, nx1 = int(np.round(shape[1]/2 - Npx)), int(np.round(shape[1]/2 + Npx)), 
ny0, ny1 = int(np.round(shape[0]/2 - Npy)), int(np.round(shape[0]/2 + Npy))

img_ = img[ny0:ny1,nx0:nx1]
extent = [-dx*Npx, dx*Npx, -dy*Npy, dy*Npy]
im = ax.imshow(img_/np.mean(img_), origin="lower", cmap="inferno", extent=extent, vmin=vmin, vmax=vmax)

H,x_edges = np.histogram(img_, bins=Tbins, range=Trange, density=True)
ax_.stairs(H*100, x_edges, label="e")
ax_.set_xlabel("Brightness temperature [K]")
ax_.set_ylabel("% of area")
ax_.legend()

ax.text(0.9,0.9, letter_list[4], color="white", size=20, transform=ax.transAxes)

ax.xaxis.set_label_position("top")
ax.xaxis.tick_top()
ax.set_xlabel("X [arcsec]")

ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_ylabel("Y [arcsec]")

#AXS = np.append(axs, ax)
#fig.colorbar(im, ax=AXS.ravel().tolist(), location="right", fraction=0.046, pad=0.02, label="Tb/mean(Tb)")
# [x position, y position, x width, y width
cax = fig.add_axes([0.98, 0.11, 0.025, 0.77])
fig.colorbar(im, cax=cax, fraction=0.046, pad=0.02, label=r"$T_B/\langle T_B \rangle$")

figname = "ALMA_comparison.pdf"
plt.savefig("figures/"+figname, bbox_inches="tight")
