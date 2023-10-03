import numpy as np
import matplotlib.pyplot as plt

import os

import h5py

input_path = "/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/"
output_path = "/mn/stornext/d19/RoCS/alma/emissa_sim/linfor3D/outhdf/"

T_list = [6500, 5770, 5000, 4000, 3200]
B_list = [0, 100, 200, 500]


### Most relevant input files:

input_files = {}
for T in T_list:
    input_files[T] = {}
    
input_files[6500][500] = input_path + "d3t65g45c_v500G_n001_it004_full_art.h5"
input_files[6500][200] = input_path + "d3t65g45c_v200G_n001_it004_full_art.h5"
input_files[6500][100] = input_path + "d3t65g45c_v100G_n001_it004_full_art.h5"
input_files[6500][0]   = input_path + "d3t65g45c_v000G_n002_it003_full_art.h5"
#input_files[6500][0]   = input_path + "d3t65g45c_v000G_end_art.h5"

input_files[5770][500] = input_path + "d3gt57g44v500_chro01_n016_it001_full_art.h5"
input_files[5770][200] = input_path + "d3gt57g44v200_chro01_n016_it001_full_art.h5"
input_files[5770][100] = input_path + "d3gt57g44v100_chro01_n016_it003_full_art.h5"
input_files[5770][0]   = input_path + "d3t57g44c_v000G_n013_it000_full_art.h5"
#input_files[5700][0]   = input_path + "d3t57g44c_v000G_n019_it000_full_art.h5"

input_files[5000][500] = input_path + "d3t50g45c_v500G_n005_it004_full_art.h5"
input_files[5000][200] = input_path + "d3t50g45c_v200G_n005_it004_full_art.h5"
input_files[5000][100] = input_path + "d3t50g45c_v100G_n002_it004_full_art.h5"
input_files[5000][0]   = input_path + "d3t50g45c_v000G_n013_it000_full_art.h5"

input_files[4000][500] = input_path + "d3t40g45c_v500G_n005_it004_full_art.h5"
input_files[4000][200] = input_path + "d3t40g45c_v200G_n005_it004_full_art.h5"
input_files[4000][100] = input_path + "d3t40g45c_v100G_n005_it004_full_art.h5"
input_files[4000][0]   = input_path + "d3t40g45c_v000G_n013_it000_full_art.h5"
#input_files[4000][0]   = input_path + "d3t40g45c_v000G_n012_end_art.h5"
#input_files[4000][0]   = input_path + "d3t40g45c_v000G_n013_it000_full_art.h5"

input_files[3200][500] = input_path + "mdw3d_t32g45_cm002_08_n048_end_art.h5"
input_files[3200][200] = input_path + "mdw3d_t32g45_cm002_07_n026_end_art.h5"
input_files[3200][100] = input_path + "mdw3d_t32g45_cm002_02_n040_end_art.h5"
input_files[3200][0]   = input_path + "mdw3d_t32g45_cm002_11_n029_it003_full_art.h5"

input_files[3200][500] = input_path + "mdw3d_t32g45_cm002_08_n044_end_art.h5"
input_files[3200][200] = input_path + "mdw3d_t32g45_cm002_07_n022_end_art.h5"
input_files[3200][100] = input_path + "mdw3d_t32g45_cm002_02_n036_end_art.h5"
input_files[3200][0]   = input_path + "mdw3d_t32g45_cm002_11_n025_end_art.h5"

#input_files[3200][500] = input_path + "mdw3d_t32g45_cm002_08_end_art.h5"
#input_files[3200][200] = input_path + "mdw3d_t32g45_cm002_07_end_art.h5"
#input_files[3200][100] = input_path + "mdw3d_t32g45_cm002_02_end_art.h5"
#input_files[3200][0]   = input_path + "mdw3d_t32g45_cm002_11_n029_it003_full_art.h5"
#input_files[3200][0]   = input_path + "d3t32g45c_v000G_n004_it000_full_art.h5"


### Most relevant output files:

output_files = {}
for Teff in input_files.keys():
    output_files[Teff] = {}
    for B in input_files[Teff].keys():
        output_files[Teff][B] = {}
        
text_list = ["05000", "00_30mm", "01_00mm", "03_00mm", "06_00mm", "09_00mm"]
val_list = [0.0005, 0.3, 1.0, 3.0, 6.0, 9.0]

def fill_output_files(Teff, B, filename0, filename1):
    """
    """
    for text, val in zip(text_list, val_list):
        filename = filename0 + text + filename1
        file_path = output_path + filename
        if os.path.exists(file_path):
            output_files[Teff][B][val] = output_path + filename


fill_output_files(6500, 500, "d3t65g45c_v500G_n001_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(6500, 200, "d3t65g45c_v200G_n001_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(6500, 100, "d3t65g45c_v100G_n001_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(6500,   0, "d3t65g45c_v000G_n002_it003_", "_mu1_00_linfor_3D_2.hdf")

fill_output_files(5770, 500, "d3gt57g44v500_chro01_n016_it001_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(5770, 200, "d3gt57g44v200_chro01_n016_it001_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(5770, 100, "d3gt57g44v100_chro01_n016_it003_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(5770,   0, "d3t57g44c_v000G_n013_it000_", "_mu1_00_linfor_3D_2.hdf")
#fill_output_files(5700, 100, "d3gt57g44v100_chro01_n016_it003_", "_mu1_00_linfor_3D_2.hdf")
#output_files[5700][0][0.0005] = output_path + "d3t57g44c_v000G_n019_it000_05000_mu1_00_linfor_3D_2.hdf"
#output_files[5700][0][1]      = output_path + "d3t57g44c_v000G_n019_it000_01mm_mu1_00_linfor_3D_2.hdf"
#output_files[5700][0][3]      = output_path + "d3t57g44c_v000G_n019_it000_03mm_mu1_00_linfor_3D_2.hdf"

fill_output_files(5000, 500, "d3t50g45c_v500G_n005_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(5000, 200, "d3t50g45c_v200G_n005_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(5000, 100, "d3t50g45c_v100G_n002_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(5000,   0, "d3t50g45c_v000G_n013_it000_", "_mu1_00_linfor_3D_2.hdf")

fill_output_files(4000, 500, "d3t40g45c_v500G_n005_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(4000, 200, "d3t40g45c_v200G_n005_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(4000, 100, "d3t40g45c_v100G_n005_it004_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(4000,   0, "d3t40g45c_v000G_n013_it000_", "_mu1_00_linfor_3D_2.hdf")

"""
fill_output_files(3200, 500, "mdw3d_t32g45_cm002_08_n048_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200, 200, "mdw3d_t32g45_cm002_07_n026_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200, 100, "mdw3d_t32g45_cm002_02_n040_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200,   0, "mdw3d_t32g45_cm002_11_n029_it003_", "_mu1_00_linfor_3D_2.hdf")

fill_output_files(3200, 500, "mdw3d_t32g45_cm002_08_n046_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200, 200, "mdw3d_t32g45_cm002_07_n024_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200, 100, "mdw3d_t32g45_cm002_02_n036_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200,   0, "mdw3d_t32g45_cm002_11_n025_it000_", "_mu1_00_linfor_3D_2.hdf")
"""

fill_output_files(3200, 500, "mdw3d_t32g45_cm002_08_n044_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200, 200, "mdw3d_t32g45_cm002_07_n022_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200, 100, "mdw3d_t32g45_cm002_02_n038_it000_", "_mu1_00_linfor_3D_2.hdf")
fill_output_files(3200,   0, "mdw3d_t32g45_cm002_11_n027_it000_", "_mu1_00_linfor_3D_2.hdf")




time_input_path = "/mn/stornext/d19/RoCS/svenwe/jonast/data/art/input/"
time_output_path = "/mn/stornext/d19/RoCS/jonasrth/ART/time_series_CF/"
time_input_path_test = "/mn/stornext/d19/RoCS/jonasrth/ART_input/time_series_backup/"

n_list = [13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19]
it_list = [0, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]
time_list = [1800.01, 1810.04, 1820.00, 1830.02, 1840.00, 1850.02, 1860.00, 1870.01, 1880.01, 1890.02, 
             1900.00, 1910.01, 1920.02, 1930.00, 1940.01, 1950.01, 1960.03, 1970.01, 1980.02]

time_input_files = []
time_output_files = []
time_input_files_test = []

for i,(n,it) in enumerate(zip(n_list, it_list)):
    time_output_files.append(time_output_path + f"output_n0{n}_it00{it}.h5")
    time_input_files.append(time_input_path + f"d3t57g44c_v000G_n0{n}_it00{it}_full_art.h5")
    time_input_files_test.append(time_input_path_test + f"d3t57g44c_v000G_n0{n}_it00{it}_full_art.h5")
