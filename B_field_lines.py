import numpy as np
import matplotlib.pyplot as plt

import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

def plot_streamlines(Teff, B0):
    """
    """
    
    iy = 100
    
    F = h5py.File(input_files[Teff][B0], "r")
    #f = h5py.File(output_files[Teff][B0][0.0005], "r")
    
    TR = np.array(F["tau"][0,:,iy,:])[:,::-1]
    iTR1 = np.argmin(abs(TR-1), axis=1)
    iTRC = np.argmin(abs(TR-1e-5), axis=1)
    i_s = round(np.mean(iTR1))

    T = np.array(F["temperature"][0,:,iy,:])[:,::-1].T
    Pg = np.array(F["Pgas"][0,:,iy,:])[:,::-1].T
    bx = np.array(F["bx"][0,:,iy,:])[:,::-1].T
    by = np.array(F["by"][0,:,iy,:])[:,::-1].T
    bz = np.array(F["bz"][0,:,iy,:])[:,::-1].T
    b = np.sqrt(bx**2 + by**2 + bz**2)
    Pb = b**2/(8*np.pi)
    beta = Pg/Pb


    z = np.array(F["z"][0,0,0,:])[::-1]*1e-8
    z = z - z[i_s]

    dx = np.array(F["dx"][0])*1e-8
    dz = abs(z[-1] - z[-2])

    z_int = np.arange(T.shape[0])
    x_int = np.arange(T.shape[-1])

    x = x_int*dx 
    z = z_int*dz + z[0]
    
    fs = np.array([10,7])*0.8
    fig,ax = plt.subplots(figsize=(fs[0],fs[1]))

    extent = [x[0],x[-1],z[0], z[-1]]
    im = ax.imshow(np.log10(b), origin="lower", extent=extent)
    plt.colorbar(im, fraction=0.0335, pad=0.04, label="$\log_{10}($Magnetic field strength [G]$)$")

    ax.streamplot(x, z, bx, bz, color="white", density=2)

    k_beta1 = np.argmin(abs(beta-1), axis=0)
    z_beta1 = z[k_beta1]

    #ax.plot(x, z_beta1, color="red", ls="dashed")
    ax.plot(x, z[k_beta1], color="red", ls="dotted")
    ax.plot(x, z[iTR1], color="red", ls="dashed")
    ax.plot(x, z[iTRC], color="red", ls="dashdot")
    
    #ax.legend()
    ax.set_xlim([extent[0], extent[1]])
    ax.set_ylim([extent[2], extent[3]])
    ax.set_xlabel("x [Mm]")
    ax.set_ylabel("z [Mm]")
    
    figname = "Bfield_fieldlines.pdf"
    plt.savefig("figures/" + figname, bbox_inches="tight")
    
plot_streamlines(5000, 200)