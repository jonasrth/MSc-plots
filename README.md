# MSc-plots

### Python code for MSc thesis figures

---

Python scripts visualise data from large 3D solar atmosphere simulations (simulated by CO5BOLD), and data from radiative transfer calculations performed on these simulations (using [ART](https://github.com/SolarAlma/ART), Linfor3D, [RH 1.5D](https://github.com/ITA-Solar/rh)). The data was provided in HDF5 files, and is not publically available. The figures produced by the python scripts can be found in the "figures" folder.

---

**Contents:**

- **figures** : Contains figures produced by pythons scripts in PDF format. The .pdf files have the same name as the scripts producing them (usually).
- **files** : Some pickle files with useful derived (secondary) data (Mean integrated optical depth, mean pressure scale height, mean brightness/gas temperature difference).
- **2DBslice.py, 2DTslice.py** : Produces figures showing vertical slices of magnetic field strength (B) and gas temperature (T), for all model atmospheres.
- **2D_Bhist.py, 2D_Thist.py** : 2D histograms of vertical magnetic field strength stratification (B) and gas temperature stratification (T)
- **ALMA_comp.py** : Comparing model data with ALMA observations.
- **B_field_lines.py** : Vertical slice of magnetic field strength including field line visualisation.
- **CF_time.py** : Evolution of contribution function over time (for a single column in the model atmosphere)
- **SED_EMISSA_X.py**: Comparing calculated spectral energy distribution (SED) from model atmospheres with observed stellar SEDs compiled in connection with [EMISSA](https://www.mn.uio.no/rocs/english/projects/emissa/index.html) project.
- 
