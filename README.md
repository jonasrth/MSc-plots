# MSc-plots

### Python code for MSc thesis figures

---

***Data:***

Python scripts visualise data from large 3D solar atmosphere simulations (simulated by CO5BOLD), and data from radiative transfer (RT) calculations performed on these simulations (using [ART](https://github.com/SolarAlma/ART), Linfor3D, [RH 1.5D](https://github.com/ITA-Solar/rh)). The data was provided in HDF5 files, and is not publically available. The figures produced by the python scripts can be found in the "figures" folder. Simulations had dimensions of (nx, ny, nz), where nx,ny,nz are the numbers of grid cells in horizontal (x,y) and vertical(nz) direction. Usually typical size around (500,500,300).

Simulation data:
- *ID {Dimension} : [unit] Description*
- z {nz}     : [cm] Height of gridcells along z-axis
- dx, dy {SCALAR} : [cm] Horisontal extent of one gridcell
- Pgas {nx,ny,nz} : [dyn/cm^2] Gas pressure
- dens {nx,ny,nz} : [g/cm^3] Mass density
- xne {nx,ny,nz} : [1/cm^3] Electron number density
- temperature {nx,ny,nz} : [K] Gas temperature
- kappa {nx,ny,nz} : [1/cm] Mean Rosseland opacity
- tau {nx, ny, nz} : [1] Mean Rosseland optical depth
- bx, by, bz {nx,ny,nz} : [G] Magnetic field strength, in directions x,y, and z
- vx, vy, vz {nx,ny,nz} : [cm/s] Bulk velocity, in directions x,y, and z

Radiative transfer data, ART:
- *ID {Dimension} : [unit] Description*
- Stokes_I {1, nx, ny, nw}: [erg cm^(-2) s^(-1) st^(-1) Hz^(-1)] Calculated intensity
- CFunc {1, nx, ny, nz} : [erg cm^(-3) s^(-1) st^(-1) Hz^(-1)] Contribution function (only calculated for one wavelength)
- Tau1 {1, nx, ny, nw} : [cm] Height of optical depth unity
- Wavelength {nw} : [angstrom] Wavelength for which the intensity was calculated

Radiative transfer data, Linfor3D:
- *ID {Dimension} : [unit] Description*
- intensity {nx, ny} : [erg cm^(-2) s^(-1) st^(-1) angstrom^(-1)] Calculated intensity
- contfunc {nx,ny,nz} : [erg cm^(-3) s^(-1) st^(-1) angstrom^(-1)] Contribution function
- contfctz {nz} : [cm] Height scale of contribution function (z = 0 at height of Rosseland optical depth unity)
- wavelength {1} : [angstrom] Wavelength for which the intensity was calculated


---

**Contents:**

- **figures** : Contains figures produced by pythons scripts in PDF format. The .pdf files have the same name as the scripts producing them (usually).
- **files** : Some pickle files with useful derived (secondary) data (Mean integrated optical depth, mean pressure scale height, mean brightness/gas temperature difference).
- **2DBslice.py, 2DTslice.py** : Produces figures showing vertical slices of magnetic field strength (B) and gas temperature (T), for all model atmospheres.
- **2D_Bhist.py, 2D_Thist.py** : 2D histograms of vertical magnetic field strength stratification (B) and gas temperature stratification (T)
- **ALMA_comp.py** : Comparing model data with ALMA observations.
- **B_field_lines.py** : Vertical slice of magnetic field strength including field line visualisation.
- **CF_time.py** : Evolution of contribution function over time (for a single column in the model atmosphere)
- **SED_EMISSA_Sun/mag/plots.py**: Comparing calculated spectral energy distribution (SED) from model atmospheres with observed stellar SEDs compiled in connection with [EMISSA](https://www.mn.uio.no/rocs/english/projects/emissa/index.html) project.
- **TR_dict.py** : Produces files/TR_dict.pkl, which contains mean optical depth as a function of height for all model atmospheres.
- **T_hist_chsp.py** : Histograms of the gas temperature distribution in the chromospheres of the model atmospheres.
- **Tb_vs_Tg.py** : Brightness temperature (Tb) vs. Gas temperature (Tg) scatter plots.
- **TbvsTg_contour.py** : Brightness temperature (Tb) vs. Gas temperature (Tg) contour plots.
- **TgTb_error.py** : Plots showing difference between Tb and Tg for different wavelengths, for all the model atmospheres.
- **Thist_chsp.py** : Old version of T_hist_chsp.py
- **Tslice.py, TSlice_mag.py** : Old versions of 2DTSlice.py and 2DBSlice.py
- **X_strat.py** : Plots stratification of variable X (density, pressure, temperature, etc.)
- **codecomp_cf.py, codecomp_map.py** : Plots comparing contribution functions (CF) or brightness temperature maps from the different RT codes.
- **int_img.py, int_imgs.py** : Plots intensity images from all the model atmospheres.
- **material.py** : Sorts filenames of all relevant data files used in MSc project into dictionaries. Imported and used by most of the other python scripts to help collect the relevant data.
- **mean_CF.py** : Plots mean contribution functions (CF) over entire horizontal extent of the simulations, as a function of height.
- **obssun_sst1.py** : Plots image of solar observation by SST.
- **plot_params.py** : Plots model atmosphere parameters in surface gravity vs. effective temperature plot, next to stellar observations from EMISSA.
- **t32_test.py** : Plots vertical slices of gas velocity magnitude for the coolest models (Effective temperature of 3200 K)
- **temp_dens_B_strat.py** : Plots temperature, density, and magnetic field strength stratification
- **time_quantity_Eklund.py** : Plots quantity (temperature, density, etc.) in a column as a function of time. Similar to plots showed in paper ([Eklund, H. 2021](https://ui.adsabs.harvard.edu/abs/2021RSPTA.37900185E/abstract)).
- **val81c.py** : Plots traditional VAL (Vernazza, Avrett, Loeser) model of solar gas temperature- and density stratification.
- **zt1_hist.py** : Histograms showing distribution of height of optical depth unity for a set of wavelengths from a single model atmosphere.
- **zt1_hist_mag.py** : Histograms showing height of optical depth unity for a single wavelength, in model atmospheres with different initial magnetic field strength.
