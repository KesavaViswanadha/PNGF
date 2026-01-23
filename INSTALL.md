# INSTALL.md

## C++ and Fortran build tools

First, install `g++`, `make`, and `gfortran`. `gfortran` is used only to compile MUMPS. For this purpose, a different C and/or Fortran compiler may need to be used, depending on your system.

## Eigen

Download [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and extract the compressed archive. 

Only the header files, at the `Eigen/` subdirectory, are necessary, and these can be placed in any convenient location, such as `~/.local/Eigen/`

## METIS

METIS is used by MUMPS for ordering. METIS may be obtained from [GitHub](https://github.com/KarypisLab/METIS), but [GKLib](https://github.com/KarypisLab/GKlib/tree/master) is required as well. From any convenient directory,

1. Extract the compressed archive (if necessary)
2. Set the source files to use double precision
    ```
    sed -i "43s/32/64/" metis-5.1.0/include/metis.h
    ```
    On MacOS, you may need to run:
    ```
    sed -i '' "43s/32/64/" metis-5.1.0/include/metis.h
    ```
3. `cd metis-5.1.0/`
4. `make config cc=gcc prefix=~/.local/metis`, where `~/.local/metis` is any convenient location. The minimum CMake version in `metis-5.1.0/CMakeLists.txt` may need to be changed to 3.5 for modern CMake to work.
5. `make install`
This will copy the METIS binaries, the header, and the library in `bin/`, `include/` and `lib/` folders to `~/.local/metis`.

On MacOS, METIS may also be installed via Homebrew.

## BLAS and LAPACK

### Linux

Download the [Intel OneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html), following the instructions for your Linux distribution. Instructions for Ubuntu using the APT package manager are below:

1. Add APT sources
```
sudo apt update
sudo apt install -y gpg-agent wget
```
2. Download the key to the system keyring
```
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB  | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
```
3. Add the signed entry to APT sources and configure to use Intel's repository:
```
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
```
4. `sudo apt update`
5. Install with `sudo apt install intel-oneapi-hpc-toolkit`

On a HPC environment, the OneAPI and related (e.g. Intel MKL) environment variables might be maintained using the [Environment Modules](https://modules.sourceforge.net/) system. In this case, your `.bashrc` should have something like the following
```
export MODULEPATH=${MODULEPATH}:/opt/intel/oneapi/modulefiles
module load mkl/2025.0
```
If not, prior to compiling MUMPS, run
```
source /opt/intel/oneapi/2025.0/oneapi-vars.sh
```
to set up the environment variables. This may also be appended to your `.bashrc` if desired.

### MacOS

The [Accelerate framework](https://developer.apple.com/documentation/accelerate) is preinstalled on ARM-based Macs (Apple Silicon). 

When compiling `pngf-opt`, ensure that the Makefile has `-framework Accelerate` in the `LDFLAGS` in order to link to it.

## MUMPS

Download [MUMPS](https://mumps-solver.org/). Uncompress and copy the contents of the archive to somewhere convenient, say `~/MUMPS_5.7.3/`.

Next, copy either `Makefile_Linux.inc` or `Makefile_Mac.inc` from the `mumps-install/` directory in this repository to `~/MUMPS_5.7.3/`, depending on whether you are working on Linux or MacOS. Rename this file to just `Makefile.inc`. Then open the file and modify the path names in the following as appropriate:
- `LMETISDIR` and `IMETIS` with the `lib/` and `include/` directory of METIS
- On Linux, the path in `RPATH_OPT` to where the MUMPS library directory is (where the shared libraries will be placed after compilation), e.g. `~/MUMPS_5.7.3/lib/`
- `MKLROOT`, `LAPACK`, and `SCALAP` if necessary for Intel MKL (Linux) or Accelerate (MacOS)
- 'LIBPAR' and 'LIBBLAS' as appropriate for Open MPI, if needed
- `OPTF`, and `OPTC` with the path to Open MPI `lib/` and `include/` directories as appropriate.
After saving `Makefile.inc`, run either `make allshared` to make shared libraries of MUMPS for single and double precision for real and complex variables. Alternatively, run `make cshared` to build only the single-precision complex variable shared library, which is the only one used im `acmefdfd`.

On MacOS, the C compiler `CC` may need to be replaced with e.g. `gcc-15`, since Apple Clang does not support `-fopenmp`. Open MPI and SCALAPACK can be installed via Homebrew.

Once MUMPS is compiled,
- If on Linux, add the `lib/` path to the `LD_LIBRARY_PATH` environment variable. For example, you may add the following to your `.bashrc.`:
```
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/MUMPS_5.7.3/lib"
```
- If on MacOS, add it to the `DYLD_LIBRARY_PATH`; for example, in your `.zshrc` or `.bashrc` (depending on which shell you are using),
```
export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$HOME/MUMPS_5.7.3/lib"
```

## acmefdfd

Edit the include and lib flags in the Makefile in the `acmefdfd/` directory with the paths for Open MPI, Eigen, MUMPS, and (if on Linux) Intel OneAPI. Then `make`.

## pngf-opt

Edit the include and lib flags in the Makefile in the `pngf-opt/` directory with the paths for Open MPI and (if on Linux) Intel OneAPI and MKL or (if on MacOS) Accelerate. Then `make`.
