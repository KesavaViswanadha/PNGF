# Precomputed Numerical Green Function Method:<br/>Near real-time full-wave electromagnetic inverse design

This repository contains code to precompute the numerical Green's function (NGF) matrix and perform optimization utilizing the PNGF method and direct binary search as described in the paper "Near real-time full-wave inverse design of electromagnetics devices" by Sun, et. al.

## Directory listing

```
/
|--- acmefdfd/ - Generation of the FDFD formulation system matrices and APF precomputation
|--- pngf-opt/ - Optimization using PNGF with the NGF matrix precomputed by APF
|--- mumps-install/ - Installation helpers for the MUMPS solver package
```

## Dependencies

A Linux (with Intel processor) or MacOS (Apple Silicon) system is assumed; testing was done by the authors using Ubuntu 24.04 and MacOS 15 Sequoia. Below is a list of all dependencies:

- C++ and Fortran build tools
- [Open MPI](https://www.open-mpi.org/)
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) linear algebra library
- [METIS](https://github.com/KarypisLab/METIS) partitioning library
- [MUMPS](https://mumps-solver.org/index.php) parallel sparse direct solver
- BLAS and LAPACK linear algebra packages
  - If on Linux, the [Intel OneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html), or
  - If on MacOS, the [Accelerate framework](https://developer.apple.com/documentation/accelerate)

For detailed installation instructions, please refer to `INSTALL.md`.

## Usage

Once the `acmefdfd` and `pngf-opt` executables have been built, 

1. `cd acmefdfd/`
2. Run the `run_acmefdfd.sh` script, which calls `acmefdfd` to generate G matrices for each of the frequencies of optimization. By default, these matrices are placed in the `pngf-opt/` directory. 
3. `cd ../pngf-opt`
4. Run the `pngf-opt` executable to perform optimization. A log file with the objective function history will be generated, and the initial and final designs are printed out.

The PNGF optimizer is configured by default for the substrate antenna design example, with a set RNG seed that generates the design reported in the paper. The executable may be recompiled to use a time-based seed to produce random designs.

For more details, please refer to the `README.md` in each of the `acmefdfd/` and `pngf-opt/` directories.

## Details

The code is configured by default to optimize a wideband 30GHz substrate antenna that uses a 1.38mm-thick dielectric (epsilon_r = 3.5) with copper clad (modeled as perfect electrical conductor material). The bottom is completely filled in as a ground plane, and the top is the optimization region. The optimization region comprises 21x21 tiles, where each tile is 0.5mm by 0.5mm and is formed by the faces of 3x3 adjoining Yee cells (that is, each tile has 4x3 = 12 electric field x-components and 3x4 = 12 electric field y-components). A x-directed lumped port is defined at the center tile; the port is made up of six x-components on a 3x1 area of Yee cell faces centered at the tile.

## Notes

- Ensure that the machine used to run `acmefdfd` to generate the numerical Green function matrices has sufficient RAM. The peak memory usage for the configured example is around 50 GB. 

## Authors

Jui-Hung Sun [1], Mohamed Elsawaf [1], Yifei Zheng [1], Ho-Chun Lin [1], Chia Wei Hsu [1], Constantine Sideris [1,2]\
[Analog/RF Integrated Circuits, Microsystems, and Electromagnetics Laboratory](https://acme.stanford.edu/)\
[Optics in Complex Systems](https://sites.usc.edu/hsugroup/)\
[1] Ming Hsieh Department of Electrical and Computer Engineering, University of Southern California\
[2] Department of Electrical Engineering, Stanford University

## License

Copyright © 2025, 2026, Constantine Sideris and Jui-Hung Sun. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 
