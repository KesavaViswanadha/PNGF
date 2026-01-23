// ACMEFDFD V1.00
//
// Simple FDFD solver for producing PNGF
// precomputation matrices.
// 
//
// Copyright (c) 2025, 2026, Constantine Sideris (sideris@stanford.edu) and Jui-Hung Sun
// (juihungs@usc.edu)
// 
// This program is free software: you can redistribute it and/or modify it under the terms 
// of the GNU Affero General Public License as published by the Free Software Foundation, 
// either version 3 of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT ANY 
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License along with 
// this program. If not, see <https://www.gnu.org/licenses/>. 
//
//
//
// Constantine Sideris (sideris@stanford.edu)
// Jui-Hung (Ray) Sun (juihungs@usc.edu)
//
// This FDFD solver is configured to produce
// the PNGF G matrices needed by the PNGF
// optimizer code to inverse design the
// Substrate Antenna (SA) in the PNGF paper.
//
// To produce all 5 G matrices for frequencies
// 25GHz, 27.6GHz, 30GHz, 32.5GHz, 35GHz, run
// acmfdfd with the following parameters:
//
// ./acmefdfd Gmat_sub_01.bin 25.0E9
// ./acmefdfd Gmat_sub_02.bin 27.5E9
// ./acmefdfd Gmat_sub_03.bin 30.0E9
// ./acmefdfd Gmat_sub_04.bin 32.5E9
// ./acmefdfd Gmat_sub_05.bin 35.0E9
//

#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <list>
#include <forward_list>
#include <functional>
#include <random>
#include <fstream>

#include <thread>
#include <chrono>
#include <time.h>

#include <stdio.h>
#include <math.h>
#include <complex>

#include <omp.h>

#include <Eigen>

#include "mpi.h"

// Use single precision float for the APF
// matrix factorization:
#define USE_MUMPS_SINGLE

// Use double precision elsewhere
typedef double precision;

#ifdef USE_MUMPS_SINGLE
#include "cmumps_c.h"
typedef float mumps_precision;
#else
#include "zmumps_c.h"
typedef double mumps_precision;
#endif

// Constants for MPI required by MUMPS
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

using namespace std;
using namespace Eigen;

// `_d` literal to cast to double
double operator"" _d(long double v)
{
    return (double)v;
}

// Convenient type aliases
using ArrayXr = std::vector<precision>;
using ArrayXcr = std::vector<std::complex<precision>>;
static const std::complex<precision> I1 = std::complex<precision>(0.0, 1.0);

// Physical constants
#define eps0            (8.8541878188E-12)
#define mu0             (1.25663706127E-6)
#define c0              (2.997924579998211E8)
#define eta0            (3.767303134118051e+02)

#include "acmetime.hpp"

// define precision prior to including these
#include "mode.hpp"
#include "apf_solver.hpp"
#include "nf2ff.hpp"
#include "matutil.hpp"
#include "acmefdfd.hpp"

int main (int argc, char **argv)
{    
    int max_threads;
    acmetime mytime;

    // Size of the optimization region in terms of
    // number of tiles. Note that each tile consists 
    // of 3x4 and 4x3 y directed Yee cells. The number of
    // Yee cells used to discretize each tile is
    // hardcoded in this proof-of-concept code for
    // the sake of simplicity, but it can be easily
    // changed to anything else. 
    const int tilemap_xsz = 21;
    const int tilemap_ysz = 21;
    // X and Y size of the dimensions in terms of the
    // number of 3dx x 3dy tiles:
    const int substrate_xsz = 21;
    const int substrate_ysz = 21;
    // This might be a bit confusing, but the height is defined
    // in terms of number of Yee cells instead.
    const int substrate_height = 25;
    const std::complex<precision> substrate_epsr = 3.5_d;
    const precision nf2ff_theta = 0.0_d;
    const precision nf2ff_phi = 0.0_d;

    std::string Gmat_fname;

    precision solve_freq = 30.0E9; // Default is 30GHz.
    precision lambda_design = 0.01; // Default is 10mm.

    // set OMP number of cores:
    max_threads = omp_get_max_threads();

    if (argc == 1)
    {
        // No parameters passed, output file is just "Gmat.bin":
        Gmat_fname = "Gmat.bin";
    } else if (argc == 2)  {
        // Only one argument. Assume that it is the filename:
        Gmat_fname = std::string(argv[1]);
    } else if (argc == 3) {
        Gmat_fname = std::string(argv[1]);
        solve_freq = atof(argv[2]);
    } else if (argc == 4) {
        Gmat_fname = std::string(argv[1]);
        solve_freq = atof(argv[2]);
        max_threads = atoi(argv[3]);
    } else if (argc > 4) {
        Gmat_fname = std::string(argv[1]);
        solve_freq = atof(argv[2]);
        max_threads = atoi(argv[3]);
        lambda_design = atof(argv[4]);
    } 
    
    precision sub_xsz = (precision)(substrate_xsz*3)*lambda_design/120.0_d+1.0E-8;
    precision sub_ysz = (precision)(substrate_ysz*3)*lambda_design/120.0_d+1.0E-8;
    precision sub_zsz = (precision)(substrate_height)*lambda_design/360.0_d+1.0E-8;
    
    std::cout << "Solve frequency: " << solve_freq << " Hz" << std::endl;
    std::cout << "Design lambda: " <<  lambda_design << " m" << std::endl;
    std::cout << "Output filename: " << Gmat_fname << std::endl;
    std::cout << "Max OMP threads: " << max_threads << std::endl;

    // Initialize 3D FDFD simulation environment:
    acmefdfd myFDFD(48+substrate_xsz*3,48+substrate_ysz*3,101,solve_freq, lambda_design);
    std::cout << "Domain size: " << myFDFD.Nx << " x " << myFDFD.Ny << " x " << myFDFD.Nz << std::endl;

    // Initialize a nonuniform grid taking the dielectric substrate dimensions into account:
    myFDFD.create_nonuniform_grid (lambda_design/60.0_d, lambda_design/60.0_d, lambda_design/180.0_d,
                                   lambda_design/20.0_d, lambda_design/20.0_d, lambda_design/20.0_d,
                                   substrate_xsz, substrate_ysz, substrate_height);

    // Build differential operators and set up PML:
    myFDFD.build_DeDh_mats_triplets();
    myFDFD.setup_pml(10,10,10,100.0,3.5);

    // Center of optimization region (cell indices):
    int dipo_cx = (myFDFD.Nx-1)/2;
    int dipo_cy = (myFDFD.Ny-2)/2;
    int dipo_cz = (myFDFD.Nz-2)/2;

    // Bottom left corner of optimization region (cell indices):
    int tile_start_x = dipo_cx-1-(tilemap_xsz-1)*0.5*3-1;
    int tile_start_y = dipo_cy-(tilemap_ysz-1)*0.5*3 - 1;
    // Bottom left corner of substrate region (cell indices):
    int sub_start_x = dipo_cx-1-(substrate_xsz-1)*0.5*3-1;
    int sub_start_y = dipo_cy-(substrate_ysz-1)*0.5*3 - 1;

    // Indices for filling in tiles with PEC
    int tile_x, tile_y;    

    // Clear PEC indicator arrays:
    for (int i = 0; i < myFDFD.numEx; i++)
        myFDFD.pec[0][i]=0;
    for (int i = 0; i < myFDFD.numEy; i++)
        myFDFD.pec[1][i]=0;
    for (int i = 0; i < myFDFD.numEz; i++)
        myFDFD.pec[2][i]=0;

    // add a Rogers dielectric substrate with epsr=3.5:
    precision xp, yp, zp;

    // Checks if a given point (cur_x, cur_y, cur_z) is inside the substrate
    auto is_inbox = [&] (precision cur_x, precision cur_y, precision cur_z)
    {
        if ((cur_x >= myFDFD.xc-sub_xsz) && (cur_x <= myFDFD.xc+sub_xsz) &&
            (cur_y >= myFDFD.yc-sub_ysz) && (cur_y <= myFDFD.yc+sub_ysz) &&
            (cur_z >= (myFDFD.zc-0*myFDFD.lambda0/360.0)-sub_zsz) && (cur_z <= (myFDFD.zc-0*myFDFD.lambda0/360.0)+sub_zsz))
            return 1;
        else
            return 0;
    };

    // Set up epsr array with effective epsilon values based on whether or not
    // they are on the substrate edges or corners. Simple sub-cell based volumetric 
    // averaging is used, which is not the most accurate approach in the general case;
    // however, in the special case where the field components on the outside boundaries
    // of the substrate are all tangential (rather than normal), the approach is equivalent
    // to the high-accuracy anisotropic permittivity averaging approach by S.G. Johnson.
    for (int k = 0; k < myFDFD.Nz; k++)
    {
        for (int j = 0; j < myFDFD.Ny; j++)
        {
            for (int i = 0; i < myFDFD.Nx; i++)
            {
                // check the x components:
                if ((j>0) && (k>0))
                {
                    xp = myFDFD.hx[i];
                    yp = myFDFD.y[j];
                    zp = myFDFD.z[k];

                    std::complex<precision> eff_epsr = 0.0_d;
                    for (int kk = -1; kk <= 1; kk+=2)
                    {
                        for (int jj = -1; jj <= 1; jj+=2)
                        {
                            for (int ii = -1; ii <= 1; ii+=2)
                            {
                                if (is_inbox(xp+(precision)ii*myFDFD.lambda0/60.0/3.0,
                                             yp+(precision)jj*myFDFD.lambda0/60.0/3.0,
                                             zp+(precision)kk*myFDFD.lambda0/180.0/3.0))
                                    eff_epsr += substrate_epsr;
                                else
                                    eff_epsr += 1.0_d;
                            }
                        }
                    }
                    eff_epsr /= 8.0_d;
                    myFDFD.epsr[0][i+(j-1)*(myFDFD.Nx)+(k-1)*myFDFD.Nx*(myFDFD.Ny-1)] = eff_epsr;
                }
                // check the y components:
                if ((i>0) && (k>0))
                {
                    xp = myFDFD.x[i];
                    yp = myFDFD.hy[j];
                    zp = myFDFD.z[k];

                    std::complex<precision> eff_epsr = 0.0_d;
                    for (int kk = -1; kk <= 1; kk+=2)
                    {
                        for (int jj = -1; jj <= 1; jj+=2)
                        {
                            for (int ii = -1; ii <= 1; ii+=2)
                            {
                                if (is_inbox(xp+(precision)ii*myFDFD.lambda0/60.0/3.0,
                                             yp+(precision)jj*myFDFD.lambda0/60.0/3.0,
                                             zp+(precision)kk*myFDFD.lambda0/180.0/3.0))
                                    eff_epsr += substrate_epsr;
                                else
                                    eff_epsr += 1.0_d;
                            }
                        }
                    }
                    eff_epsr /= 8.0_d;
                    myFDFD.epsr[1][(i-1)+(j)*(myFDFD.Nx-1)+(k-1)*(myFDFD.Nx-1)*myFDFD.Ny] = eff_epsr;
                }
                // check the z components:
                if ((i>0) && (j>0))
                {
                    xp = myFDFD.x[i];
                    yp = myFDFD.y[j];
                    zp = myFDFD.hz[k];

                    std::complex<precision> eff_epsr = 0.0_d;
                    for (int kk = -1; kk <= 1; kk+=2)
                    {
                        for (int jj = -1; jj <= 1; jj+=2)
                        {
                            for (int ii = -1; ii <= 1; ii+=2)
                            {
                                if (is_inbox(xp+(precision)ii*myFDFD.lambda0/60.0/3.0,
                                             yp+(precision)jj*myFDFD.lambda0/60.0/3.0,
                                             zp+(precision)kk*myFDFD.lambda0/180.0/3.0))
                                    eff_epsr += substrate_epsr;
                                else
                                    eff_epsr += 1.0_d;
                            }
                        }
                    }
                    eff_epsr /= 8.0_d;
                    myFDFD.epsr[2][(i-1)+(j-1)*(myFDFD.Nx-1)+k*(myFDFD.Nx-1)*(myFDFD.Ny-1)] = eff_epsr;
                }
            }
        }
    }

    // Fill in the metal ground plane on the bottom side of the
    // dielectric substrate.
    int pec_ground_z = 38-1; // Z index slice of the ground plane.
    int tile_map_z = 63-1;   // Z index slice of the top surface optimization region.
    for (int j = 0; j < substrate_ysz; j++)
    {
        for (int i = 0; i < substrate_xsz; i++)
        {
            tile_x = sub_start_x + i*3;
            tile_y = sub_start_y + j*3;

            for (int jj = 0; jj < 4; jj++)
            {
                for (int ii = 0; ii < 4; ii++)
                {
                    if (ii > 0)
                    {
                        myFDFD.epsr[0][(tile_x+ii)+(tile_y+jj)*(myFDFD.Nx)+pec_ground_z*myFDFD.Nx*(myFDFD.Ny-1)]=1.0-I1*1.0E12;
                    }
                    if (jj > 0)
                    {
                        myFDFD.epsr[1][(tile_x+ii)+(tile_y+jj)*(myFDFD.Nx-1)+pec_ground_z*(myFDFD.Nx-1)*myFDFD.Ny]=1.0-I1*1.0E12;
                    }
                }
            }    
        }
    }

    // compute epsr_inv:
    myFDFD.prepare_epsr_inv();

    // we are using complex epsr (see above) to represent metals inside of these PEC indicator arrays:
    std::fill(myFDFD.pec[0].begin(), myFDFD.pec[0].end(), 0);
    std::fill(myFDFD.pec[1].begin(), myFDFD.pec[1].end(), 0);
    std::fill(myFDFD.pec[2].begin(), myFDFD.pec[2].end(), 0);

    // construct the system matrix:
    myFDFD.build_system_matrix();

    std::cout << "domain center: " << myFDFD.xc << ", " << myFDFD.yc << ", " << myFDFD.zc << std::endl;

    // Build the E to EH extraction matrix since our FDFD formulation is E-field only and the NF2FF
    // transformation requires both E and H components:
    myFDFD.setup_extractEH_matrix ();

    dipo_cx = (myFDFD.Nx-1)/2;
    dipo_cy = (myFDFD.Ny-2)/2;

    // Use the PEC arrays to indicate which Yee cells on the top surface of the dielectric substrate
    // are inside the optimization tilemap grid. The  This is a bit messy, but saves using extra memory
    // since the PEC arrays of the FDFD object are not used after having constructed the system matrix operator. 
    // This works by filling out every tile inside the optimization grid, by setting the corresponding Yee cell
    // index in the PEC indicator array to 1. Thus, the PEC x and y arrays can afterwards be used to generate
    // the B and C mapping matrices, which, as described in the paper, map Jopt to Jsim and Esim to Eopt respectively.
    for (int j = 0; j < tilemap_ysz; j++)
    {
        for (int i = 0; i < tilemap_xsz; i++)
        {
            tile_x = tile_start_x + i*3;
            tile_y = tile_start_y + j*3;        
            for (int jj = 0; jj < 4; jj++)
            {
                for (int ii = 0; ii < 4; ii++)
                {
                    if (ii > 0)
                    {
                        myFDFD.pec[0][(tile_x+ii)+(tile_y+jj)*(myFDFD.Nx)+tile_map_z*myFDFD.Nx*(myFDFD.Ny-1)]=1;
                    }
                    if (jj > 0)
                    {
                        myFDFD.pec[1][(tile_x+ii)+(tile_y+jj)*(myFDFD.Nx-1)+tile_map_z*(myFDFD.Nx-1)*myFDFD.Ny]=1;
                    }
                }  
            }   
        }    
    }

    // Setup for APF:
    // B matrix is APF right hand side (1 column per optimization region point,
    // plus zero padding to make size of B^T = C)
    // C matrix is APF left hand side, which includes B^T and the projection 
    // matrix needed to compute RCS. Note that this is not the PNGF matrix
    // (I - P) + PG used in optimization.
    // Schur complement is C * A^(-1) * B
    // Note that C does not just include B^T (needed to obtain the G matrix), but also includes
    // Gobj, which in this example consists of two additional rows that are required for computing
    // the directivity of the design in the far-field at theta,phi = (0,0). Thus, 
    // rows(C)=rows(G)+rows(Gobj), where the number of rows of
    // G correspond to the number of total x and y Yee cells in the optimization region and Gobj has
    // 2 rows. Since MUMPS requires a square matrix system, the B matrix is padded with 2 additional zero
    // columns to make the overall APF matrix square.
    Eigen::SparseMatrix<std::complex<precision>> spB(myFDFD.Amat.rows(), (tilemap_xsz*3)*(tilemap_ysz*3+1) + 
                                                                         (tilemap_ysz*3)*(tilemap_xsz*3+1)+2);
    Eigen::SparseMatrix<std::complex<precision>> spC((tilemap_xsz*3)*(tilemap_ysz*3+1) + 
                                                     (tilemap_ysz*3)*(tilemap_xsz*3+1)+2, myFDFD.Amat.rows());

    // Set the size of the schur matrix block (corresponding to -C*A^(-1)*B):
    int schur_size = spB.cols();

    Eigen::SparseMatrix<std::complex<precision>> spD(schur_size,schur_size);
    Eigen::SparseMatrix<std::complex<precision>> spAPF;
        
    spB.reserve((tilemap_xsz*3)*(tilemap_ysz*3+1)+(tilemap_ysz*3)*(tilemap_xsz*3+1));
    spC.reserve((tilemap_xsz*3)*(tilemap_ysz*3+1)+(tilemap_ysz*3)*(tilemap_xsz*3+1)+myFDFD.Amat.rows()*2);
    spB.setZero();
    spC.setZero();

    // Build B and C matrices for the C A^-1 B APF computation:
    // B maps from Jopt to Jsim, C maps from Esim to Eopt.
    int Bidx=0;
    // Take care of the x Yee cell components first:
    for (int i = 0; i < myFDFD.pec[0].size(); i++)
    {
        if (myFDFD.pec[0][i])
        {
            // Due to the FDFD formulation used, the J design on the RHS must be
            // multiplied by epsr^-1/(i*eps0*w0):
            spB.insert(i, Bidx) = -myFDFD.epsr_inv.diagonal()[i]/(I1*eps0*myFDFD.w0);
            spC.insert(Bidx, i) = 1.0;
            Bidx++;
        }
    }

    // and now handle the y Yee cell components:
    for (int i = 0; i < myFDFD.pec[1].size(); i++)
    {
        if (myFDFD.pec[1][i])
        {
            spB.insert(myFDFD.numEx+i, Bidx) = -myFDFD.epsr_inv.diagonal()[myFDFD.numEx+i]/(I1*eps0*myFDFD.w0);
            spC.insert(Bidx,myFDFD.numEx+i) = 1.0;
            Bidx++;
        }
    }

    // A Near-Field to Far-Field (NF2FF) transformation is used to compute the directivity of the system at a desired
    // angle in the far-field. Note, that the NF2FF operator requires both the E and H fields in the domain on a uniformly
    // spaced grid, so the underlying mechanics of this is rather complicated, despite the fact that the final result is
    // a 2xNsim Gobj matrix. First, a matrix operator (EHfEmat) is requested from the FDFD object to extract E and H when
    // multiplying the E fields, since we use an E-only FDFD formulation to reduce the size of the FDFD system matrix A. 
    // Next, we call "generate_NFbox_matrix" of the FDFD object, which produces an interpolation matrix (NFbox_mat) that 
    // extracts E and H fields collocated at the same uniformly spaced nodes on the surface of a box enclosing the whole
    // substrate. Finally, we use the NFtoFF helper class to produce an NF2FF transform matrix (RCS_mat) that projects
    // these uniform surface E/H fields to a desired point on the infinite far-field sphere (in this example, we choose
    // the broadside direction, which corresponds to theta,phi = 0,0).
    Eigen::SparseMatrix<std::complex<precision>> NFbox_mat, RCS_mat, RCSNFbox_mat;

    int nfbox_xsz = ((substrate_xsz-1)*0.5+2)*2*3+1;
    int nfbox_ysz = ((substrate_ysz-1)*0.5+2)*2*3+1;
    int nfbox_zsz = 73;
    std::cout << "NF box size: (" << nfbox_xsz << ", " << nfbox_ysz << ", " << nfbox_zsz << ")" << std::endl;

    NFtoFF myFF(nf2ff_theta,nf2ff_theta,nf2ff_phi,nf2ff_phi,1,1,
                myFDFD.lambda0/60.0,myFDFD.lambda0/60.0,myFDFD.lambda0/60.0,nfbox_xsz,nfbox_ysz,nfbox_zsz,myFDFD.k0);
    Eigen::VectorXd rcs(1);

    myFDFD.generate_NFbox_matrix (myFDFD.xc-((precision)(substrate_xsz-1)*0.5+2.0)*myFDFD.lambda0/20.0, myFDFD.xc+((precision)(substrate_xsz-1)*0.5+2.0)*myFDFD.lambda0/20.0,
                                  myFDFD.yc-((precision)(substrate_ysz-1)*0.5+2.0)*myFDFD.lambda0/20.0, myFDFD.yc+((precision)(substrate_ysz-1)*0.5+2.0)*myFDFD.lambda0/20.0,
                                  myFDFD.zc-12.0*myFDFD.lambda0/20.0, myFDFD.zc+12.0*myFDFD.lambda0/20.0,
                                  nfbox_xsz, nfbox_ysz, nfbox_zsz, NFbox_mat);
    myFF.generate_RCS_matrix(RCS_mat);
    std::cout << "NFbox dims: " << NFbox_mat.rows() << " x " << NFbox_mat.cols() << std::endl;
    std::cout << "RCS dims: " << RCS_mat.rows() << " x " << RCS_mat.cols() << std::endl;
    RCSNFbox_mat = RCS_mat*NFbox_mat*myFDFD.EHfEmat;
    std::cout << "final dims: " << RCSNFbox_mat.rows() << " x " << RCSNFbox_mat.cols() << std::endl;

    // Insert the RCS extraction points into the left hand side matrix
    for (int i = 0; i < myFDFD.Amat.rows(); i++)
    {
        spC.insert(Bidx,i) = RCSNFbox_mat.coeff(0,i);
        spC.insert(Bidx+1,i) = RCSNFbox_mat.coeff(1,i);
    }

    spB.makeCompressed();
    spC.makeCompressed();
    spD.setZero();

    // For debugging purposes, B and C can be saved to CSC-format binary files:
    // saveSparseMatrixBinary (spB, "Bmat_sub.bin");
    // saveSparseMatrixBinary (spC, "Cmat_sub.bin");
    
    for (int i = 0; i < myFDFD.pec[0].size(); i++)
        myFDFD.pec[0][i] = 0;
    for (int i = 0; i < myFDFD.pec[1].size(); i++)
        myFDFD.pec[1][i] = 0;

    std::cout << "FDFD matrix size: " << myFDFD.Amat.rows() << "x" << myFDFD.Amat.cols() << " (" << myFDFD.numE << ")" << std::endl;

    APF_solver myAPF;

    // Schur complement
    MatrixXcd schur_out(schur_size, schur_size);

    // Build augmented APF matrix
    Eigen::SparseMatrix<std::complex<precision>> h1 = horizontalConcat(myFDFD.Amat,spB);
    Eigen::SparseMatrix<std::complex<precision>> h2 = horizontalConcat(spC,spD);
    spAPF = verticalConcat(h1,h2);

    // Run APF
    std::cout << "Running APF solver... " << std::endl;
    myAPF.provide_matrix(spAPF);
    myAPF.set_num_threads(max_threads);
    myAPF.use_L0_threads(1);

    mytime.tik();
    myAPF.solve_apf (schur_out, schur_size);
    double diff = mytime.tok();
    
    std::cout << std::string(80, '-') << std::endl;
    std::cout << "APF solve time: " << diff << " sec" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    if (!myAPF.error)
    {
        // save the G matrix (with Gobj) to a sparse matrix binary file
        FILE *gfil = fopen (Gmat_fname.c_str(), "wb");
        fwrite (schur_out.data(), 1, sizeof(std::complex<precision>)*schur_size*schur_size, gfil);
        fclose (gfil);
        return 0;
    } else {
        std::cout << "An error has occured, please check error code returned by MUMPS." << std::endl;
    }

    return 0;
}