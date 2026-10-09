# Contents: /include_misc

This folder contains auxiliary tools and libraries used in various parts of the MRI sequence generation, simulation and reconstruction pipeline.

## Subfolders

- **src_colormaps**  
  Additional colormaps for visualization of relaxation time maps (adapted from https://doi.org/10.5281/zenodo.8268884 and https://doi.org/10.5281/zenodo.1243862).

- **src_eigen**  
  Header-only linear algebra library [Eigen v3.4.0](https://eigen.tuxfamily.org).  
  Used for MEX-based EPG and Bloch simulations.

- **src_ep5**  
  Tools for MRI data visualization and reconstruction (e.g., `kspace2image`, `openadapt`, relaxation mapping, segmentation of LV ROIs).

- **src_espirit**  
  ESPIRiT implementation for coil sensitivity estimation.  
  [Source](https://people.eecs.berkeley.edu/~mlustig/Software.html),  
  [Reference: Uecker et al., MRM 2013, doi:10.1002/mrm.24751](https://doi.org/10.1002/mrm.24751)

- **src_ismrmrd**
  MATLAB implementation of the [ISMRM Raw Data Format (ISMRMRD)](https://github.com/ismrmrd/ismrmrd). This folder contains the MATLAB components required for reading ISMRMRD HDF5 datasets in OpenMRF.  
  [GitHub Repository](https://github.com/ismrmrd/ismrmrd),  
  [Reference: Inati et al., MRM 2017, doi:10.1002/mrm.26089](https://doi.org/10.1002/mrm.26089)
  
- **src_mapVBVD**  
  Functions for reading Siemens rawdata (TWIX format).  
  [GitHub Repository](https://github.com/pehses/mapVBVD) by Philipp Ehses.

- **src_mri_phantom**  
  Code for simulating structured numerical phantoms and MRI signals:
  - `basic_phantoms/`: Pre-simulated `.mat` files
  - `src_epfl/`: Realistic coil sensitivity maps ([IEEE TMI 2012](https://doi.org/10.1109/TMI.2011.2174158))
  - `src_mgram/`: Tools for signal generation and custom phantom creation

- **src_nist_values**  
  Reference values from the [NIST MRI phantom](https://doi.org/10.1002/mrm.28779).

- **src_safe**  
  Prediction of PNS in Siemens MRI systems by Filip Szczepankiewicz and Thomas Witzel: [SAFE model](https://github.com/filip-szczepankiewicz/safe_pns_prediction).

- **src_tOptGrad**  
  Calculation of time optimized gradient waveforms for arbitrary k-space trajectories by Sana Vaziri Michael Lustig: [Time Optimal Gradient Design](https://people.eecs.berkeley.edu/~mlustig/Software.html).

- **src_vds**  
  Variable-Density Spiral Design Functions by Brian Hargreaves: [VDS toolbox](http://mrsrl.stanford.edu/~brian/vdspiral/).

_Maximilian Gram: 25.08.2026_