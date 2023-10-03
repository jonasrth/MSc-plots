import numpy as np
import matplotlib.pyplot as plt

#from scipy.interpolate import RegularGridInterpolator
#from matplotlib.ticker import FormatStrFormatter

import pickle
import h5py

import sys
sys.path.append("/uio/hume/student-u19/jonasrth/Python/")

from material import input_files, output_files

#HP_dict = np.load("../../Python/files/HP_dict.pkl", allow_pickle=True)
#TR_dict = np.load("../../Python/files/TR_dict.pkl", allow_pickle=True)

TR_dict = {}

Teff_list = [6500, 5770, 5000, 4000, 3200]
B0_list = [0,100,200,500]

for Teff in Teff_list:
    TR_dict[Teff] = {}
    for B0 in B0_list:
        F = h5py.File(input_files[Teff][B0], "r")
        
        logTR = np.log10(np.mean(np.array(F["tau"][0,...]), axis=(0,1)))
        TR_dict[Teff][B0] = logTR
        
pickle.dump(TR_dict, open("files/TR_dict.pkl", "wb"))