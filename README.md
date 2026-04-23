![OpenMRF](OpenMRF_banner.png)

# Documentation

Please visit **👉[openmrf.org](https://openmrf.org)** for the official documentation & blog.

# Introduction

OpenMRF is a modular and vendor-neutral framework for Magnetic Resonance Fingerprinting (MRF) built on the open-source [Pulseq](https://pulseq.github.io) standard. It is built upon the MATLAB version of Pulseq by Kelvin J. Layton and Maxim Zaitsev ([doi:10.1002/mrm.26235](https://doi.org/10.1002/mrm.26235)). OpenMRF unifies all core components of the MRF workflow within a single MATLAB-based toolbox: flexible sequence generation, automated Bloch-based dictionary simulation, and low-rank image reconstruction. The provided tools support a wide range of contrast preparations and readouts (e.g., spiral, radial, rosette) and include integrated solutions for trajectory calibration, spin-lock modeling, slice profile simulation, and metadata storage. Designed for reproducibility and portability, OpenMRF enables easy deployment of MRF protocols across multiple scanner platforms, including Siemens, GE, Philips and United Imaging systems.

# Contents

- `include_cwru/`: Contains MRF-specific source code governed by a separate **End User License Agreement ([EULA](./include_cwru/EULA.pdf))** provided by Case Western Reserve University.

- `include_miitt/`: Contains the low-rank reconstruction code provided by the MIITT group and Jeffrey Fessler's [MIRT toolbox](https://web.eecs.umich.edu/~fessler/code/). Includes an installation script. **Do not** add this folder manually to your MATLAB path; use the `install_OpenMRF.m` script.
- `include_misc/`: Miscellaneous utilities and helper functions.

- `include_pre_sim_library/`: Library containing stored pre-simulated slice profiles, adiabatic efficiencies and compressed dictionaries.

- `include_pulseq_master/`: Copy of the official Pulseq repository ([GitHub link](https://github.com/pulseq/pulseq), v1.5.1, 09.03.2026).
- `include_pulseq_toolbox/`: Contains standard imaging readouts (cartesian, radial, spiral, rosette) combined with various preparation modules (inversion recovery, saturation, MLEV-T2, spin-lock, adiabatic spin-lock, CEST). Also includes simulation tools which can be used for dictionary generation.

- `main_sequences/`: Example Pulseq sequences and reconstruction scripts.

- `user_specifications/`: User specific definitions (automatically generated via `install_OpenMRF.m`) and MRI system specifications (create a `.csv` file for your system's gradient limits and timings).

# MRF IP Notice

Magnetic Resonance Fingerprinting (MRF) technology implemented in parts of the directory `include_cwru` is protected intellectual property owned by **Case Western Reserve University (CWRU)**. The underlying methods and related technology are subject to patent protection.

MRF is protected by the following foundational patents:

- **United States:** 8,723,518  
- **United States:** 10,241,174  
- **United States:** 10,416,259  
- **United States:** 10,627,468  
- **Europe:** EP 2897523  
- **Japan:** 6557710  
- **South Korea:** 10-1674848  
- **China:** ZL2013800592667 

The foundational scientific publication describing this technology is:

Ma D, Gulani V, Seiberlich N, Liu K, Sunshine JL, Duerk JL, Griswold MA.  
**Magnetic resonance fingerprinting.**  
*Nature.* 2013 Mar 14;495(7440):187–192.  
doi: https://doi.org/10.1038/nature11971

Use of the MRF-specific source code and associated software components contained in `include_cwru` is governed by a separate **End User License Agreement ([EULA](./include_cwru/EULA.pdf))** (`include_cwru/EULA.pdf`). The MRF technology remains the property of Case Western Reserve University and is provided subject to the terms and conditions defined in the applicable EULA.

By accessing or using the MRF functionality provided in this repository, users acknowledge and agree to comply with the terms described in the applicable EULA.

# Licenses

This repository includes third-party software distributed under their respective licenses. Please consult the [NOTICE](./NOTICE) file before use. Not all code in this repository is MIT-licensed. Users intending commercial use must review the [NOTICE](./NOTICE) file and seek appropriate permissions from the original authors.

_The OpenMRF Team: 22.03.2026_
