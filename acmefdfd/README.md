# acmefdfd

This C++ project generates the FDFD system matrices and computes the G (Schur complement) and Gobj matrices for a given design environment and optimization region configuration.

## Usage

The `acmefdfd` executable takes up to four optional arguments
1. the name of the binary file in which to store the output matrix (default: `Gmat.bin`)
2. the solution frequency (default: 30GHz)
3. the maximum number of OpenMP threads to use (default: all available threads)
4. the nominal wavelength in meters, used to define the geometry (default: 0.01)

To generate the 5 matrices needed for optimization of the substrate antenna example, which is configured by default, the `run_acmefdfd.sh` script may be used. The script calls `acmefdfd` sequentially 5 times, naming and placing the output matrix files in `pngf-opt/` for the optimizer to use. The script does not take arguments.

The output matrix is the vertical concatenation of the numerical Green function matrix G with Gobj, as described in the paper. In this particular example, Gobj is a wide matrix that has two rows, which extracts the field components required to compute the directivity at the far-field in the optimization.

## Notes

- Make sure that the machine used to run `acmefdfd` to generate the precomputation matrices has sufficient RAM. The peak memory usage for the configured example is around 50 GB. 

## Details

A high-level overview of the system construction is as follows:

1. The geometry of the optimization region and the substrate are defined. For the substrate antenna example, the optimization region is 21 by 21 tiles (as well as the substrate); and, thus,
```
const int tilemap_xsz = 21;
const int tilemap_ysz = 21;
const int substrate_xsz = 21;
const int substrate_ysz = 21;
```

2. A nonuniform finite-difference grid is created. The grid is at its finest over the substrate volume, and the Yee cell size gradually increases away from the center of the domain.

3. The discretized Maxwell difference operators DE and DH matrices are constructed.

4. Perfectly Matched Layer (PML) boundaries are defined near the bounds of the domain, for absorbing boundaries.

5. The substrate and ground planes are filled in.

6. A diagonal matrix of epsilon_r values at all points in the domain is set up and inverted. This is then used with the difference operators to create the system matrix, denoted A in the paper.

7. The APF right-hand side matrix B, which is the projection matrix that maps current excitation indices in the optimization region to electric field Yee grid points in the simulation environment, is constructed.

8. The APF left-hand side matrix, which is the vertical concatenation of B^T with a wide matrix that extracts field points outside the optimization region needed to compute the objective function in optimization (in this case, directivity), is constructed.

9. The Schur complement is computed with MUMPS and saved.
