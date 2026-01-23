// FDFD system generation class
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

#pragma once

#include <vector>
#include <iostream>
#include <complex>
#include <Eigen>

// Numeric enums for field components
#define COMPONENT_EX        0
#define COMPONENT_EY        1
#define COMPONENT_EZ        2
#define COMPONENT_HX        3
#define COMPONENT_HY        4
#define COMPONENT_HZ        5

class acmefdfd
{
// Macros for converting 3D grid indices to 1D array index
#define Ex_idx(i,j,k)           ((i)+((j)*Nx)+((k)*Nx*(Ny-1)))
#define Ey_idx(i,j,k)           (numEx+(i)+((j)*(Nx-1))+((k)*(Nx-1)*Ny))
#define Ez_idx(i,j,k)           (numEx+numEy+(i)+((j)*(Nx-1))+((k)*(Nx-1)*(Ny-1)))
#define Hx_idx(i,j,k)           ((i)+((j)*(Nx-1))+((k)*(Nx-1)*Ny))
#define Hy_idx(i,j,k)           (numHx+(i)+((j)*(Nx))+((k)*(Nx)*(Ny-1)))
#define Hz_idx(i,j,k)           (numHx+numHy+(i)+((j)*(Nx))+((k)*(Nx)*(Ny)))

#define H_idx(i,j,k)            ((i)+((j)*Nx)+((k)*Nx*Ny))

    public:
        unsigned int Nx, Ny, Nz;                    // Number of Yee cells
        std::vector<precision> x, y, z;             // grid coordinates
        std::vector<precision> hx, hy, hz;          // grid cell centers

        std::vector<precision> dxa, dya, dza;       // grid cell sizes
        std::vector<precision> dhxa, dhya, dhza;    // distance between centers

        // domain center:
        precision xc, yc, zc;

        precision dx, dy, dz;               // step size of uniform grid
        precision lambda0;                  // wavelength
        precision w0, k0;                   // angular frequency and wavenumber
        std::array<ArrayXcr,3> epsr, mur;   // permittivity/permeability in xyz
        std::array<std::vector<int>,3> pec; // is this field component PEC?
        Eigen::DiagonalMatrix<std::complex<precision>,Eigen::Dynamic> epsr_inv;
        Eigen::DiagonalMatrix<std::complex<precision>,Eigen::Dynamic> mur_inv;
        std::vector<precision> Ex, Ey, Ez;  // fields
        std::vector<precision> Hx, Hy, Hz;
        unsigned int numEx, numEy, numEz;   // E field components
        unsigned int numHx, numHy, numHz;   // H field components
        unsigned int numPMLX, numPMLY, numPMLZ;     // PML components
        unsigned int numE, numH;            // numEx + numEy + numEz

        // PML sigma arrays:
        std::vector<std::complex<precision>> sigEx, sigEy, sigEz;
        std::vector<std::complex<precision>> sigHx, sigHy, sigHz;

        // Maxwell differential operators:
        Eigen::SparseMatrix<std::complex<precision>> DE, DH;

        // FDFD system matrix:
        Eigen::SparseMatrix<std::complex<precision>> Amat;
        // Multiply E fields by this matrix to extract the H fields:
        Eigen::SparseMatrix<std::complex<precision>> EHfEmat;
        
        acmefdfd (int mNx, int mNy, int mNz, precision freq, precision lambda_design)
        {
            Nx = mNx;
            Ny = mNy;
            Nz = mNz;
            
            // field components do not include those on the PML boundaries
            numEx = Nx*(Ny-1)*(Nz-1);
            numEy = (Nx-1)*Ny*(Nz-1);
            numEz = (Nx-1)*(Ny-1)*Nz;
            numHx = (Nx-1)*Ny*Nz;
            numHy = Nx*(Ny-1)*Nz;
            numHz = Nx*Ny*(Nz-1);

            // total number of components
            numE = numEx + numEy + numEz;
            numH = numHx + numHy + numHz;

            w0 = 2*M_PI*freq;
            k0 = w0/c0; //(2.0_d*M_PI)/lambda0;

            // For a multi-frequency optimization problem, one PNGF matrix
            // is required for each frequency being considered. Each of these
            // PNGF matrices must represent the same physical design
            // structure despite being discretized at a different w0/k0, so
            // keep lambda0 fixed to the design point dimensions. In our
            // susbtrate antenna design example, we choose lambda_design=10mm,
            // which corresponds to a center design frequency of approx. 30GHz.
            // All of the other length-based units (i.e., dx, dy, dz) are derived
            // based on this design lambda:
            lambda0 = lambda_design;
            std::cout << "design lambda: " << lambda0 << std::endl;
    
            dx=lambda0/20.0_d;
            dy=lambda0/20.0_d;
            dz=lambda0/20.0_d;

            x.resize(Nx+1);
            y.resize(Ny+1);
            z.resize(Nz+1);
            hx.resize(Nx);
            hy.resize(Ny);
            hz.resize(Nz);
            dxa.resize(Nx);
            dya.resize(Ny);
            dza.resize(Nz);
            dhxa.resize(Nx-1);
            dhya.resize(Ny-1);
            dhza.resize(Nz-1);

            epsr[0].resize(numEx);
            epsr[1].resize(numEy);
            epsr[2].resize(numEz);
            mur[0].resize(numHx);
            mur[1].resize(numHy);
            mur[2].resize(numHz);

            pec[0].resize(numEx);
            pec[1].resize(numEy);
            pec[2].resize(numEz);

            epsr_inv.resize(numE);
            mur_inv.resize(numH);

            // initialize permittivity and PEC indicators
            for (int i = 0; i < numEx; i++)
            {
                epsr[0][i] = 1.0;
                pec[0][i] = 0;
            }
            for (int i = 0; i < numEy; i++)
            {
                epsr[1][i] = 1.0;
                pec[1][i] = 0;
            }
            for (int i = 0; i < numEz; i++)
            {
                epsr[2][i] = 1.0;                
                pec[2][i] = 0;
            }

            // initialize permeability
            for (int i = 0; i < numHx; i++)
                mur[0][i] = 1.0;
            for (int i = 0; i < numHy; i++)
                mur[1][i] = 1.0;
            for (int i = 0; i < numHz; i++)
                mur[2][i] = 1.0;             

           std::cout << "FDFD class constructed." << std::endl;
        }

        // substrate_xsz and substrate_ysz are passed as parameters because they are used to determine how to produce
        // the nonuniform grid. The grid dx/dy/dz are setup to be dx_min x dy_min x dz_min inside the substrate, whereas
        // they are grown to dx_max x dy_max x dz_max in the airbox outside.
        void create_nonuniform_grid (precision dx_min, precision dy_min, precision dz_min,
                                     precision dx_max, precision dy_max, precision dz_max,
                                     int substrate_xsz, int substrate_ysz, int substrate_height)
        {

            int nu_boundary_x = (substrate_xsz*3+3)/2;
            int nu_boundary_y = (substrate_ysz*3+3)/2;
            int nu_boundary_z = (substrate_height+3)/2;

            std::cout << "NU boundaries: " << nu_boundary_x << " " << nu_boundary_y << " " << nu_boundary_z << std::endl;
            
            // create a nonuniform grid
            // initial finest grid spacing at the center of the grid, where
            // the optimization region is, and grow spacing further from center
            dx = dx_min;
            dxa[(Nx-1)/2] = dx;
            for (int i = 1; i <= (Nx-1)/2; i++)
            {
                if (i>nu_boundary_x)
                {
                    // start growing dx:
                    dx = dx*1.2;
                    if (dx>(dx_max))
                        dx = dx_max;
                }
                dxa[(Nx-1)/2+i] = dx;
                dxa[(Nx-1)/2-i] = dx;
            }
            dy = dy_min;
            dya[(Ny-1)/2] = dy;
            for (int i = 1; i <= (Ny-1)/2; i++)
            {
                if (i>nu_boundary_y)
                {
                    // start growing dx:
                    dy = dy*1.2;
                    if (dy>(dy_max))
                        dy = dy_max;
                }
                dya[(Ny-1)/2+i] = dy;
                dya[(Ny-1)/2-i] = dy;
            }
            dz = dz_min;
            dza[(Nz-1)/2] = dz;
            for (int i = 1; i <= (Nz-1)/2; i++)
            {
                if (i>nu_boundary_z)
                {
                    // start growing dx:
                    dz = dz*1.2;
                    if (dz>(dz_max))
                        dz = dz_max;
                }
                dza[(Nz-1)/2+i] = dz;
                dza[(Nz-1)/2-i] = dz;
            }

            // update the grid coordinates:
            x[0] = 0;
            for (int i = 1; i < Nx+1; i++)
            {
                x[i] = x[i-1] + dxa[i-1];
                hx[i-1] = 0.5_d*(x[i]+x[i-1]);
            }
            for (int i = 0; i < Nx-1; i++)
                dhxa[i] = hx[i+1]-hx[i];

            y[0] = 0;
            for (int i = 1; i < Ny+1; i++)
            {
                y[i] = y[i-1] + dya[i-1];
                hy[i-1] = 0.5_d*(y[i]+y[i-1]);
            }
            for (int i = 0; i < Ny-1; i++)
                dhya[i] = hy[i+1]-hy[i];

            z[0] = 0;
            for (int i = 1; i < Nz+1; i++)
            {
                z[i] = z[i-1] + dza[i-1];
                hz[i-1] = 0.5_d*(z[i]+z[i-1]);
            }
            for (int i = 0; i < Nz-1; i++)
                dhza[i] = hz[i+1]-hz[i];

            // find domain center locations
            xc = (x[0]+x[Nx])*0.5_d;
            yc = (y[0]+y[Ny])*0.5_d;
            zc = (z[0]+z[Nz])*0.5_d;
        }

        // helper function for interpolator:
        // interpolant should be: f[idx]*a + f[idx+1]*(1-a)
        // returns: 
        //    0 - x is out of range
        //    1 - x was found exactly within tol of an index in xvals
        //    2 - x is in between two indicies of xvals
        int find_index (std::vector<precision> &xvals, precision x, precision *ratio, int *idx)
        {
            const precision tol=1.0E-15;
            // loop through the whole array of Cartesian coordinates
            // to find in between which two the target point x lies:
            for (int i = 0; i < xvals.size()-1; i++)
            {
                // check if x is equal to an existing element within tolerance tol:
                if (fabs(xvals[i]-x)<tol)
                {
                    (*idx) = i;
                    (*ratio) = 1.0_d;
                    return 1;
                // check the same thing for the other end of the interval:
                } else if (fabs(xvals[i+1]-x)<tol) {
                    (*idx) = i+1;
                    (*ratio) = 1.0_d;
                    return 1;
                // check if x is inbetween xvals[i] and xvals[i+1]:
                } else if ( ((x-xvals[i])>0) && ((x-xvals[i+1])<0) ) {
                    (*idx) = i;
                    (*ratio) = (xvals[i+1]-x)/(xvals[i+1]-xvals[i]);
                    return 2;
                }
            }
            return 0;
        }

        // Generate an interpolation matrix that multiplies against the full FDFD (E/H) fields
        // vector and produces the tangential E/H fields collocated over the same grid points
        // on the surfaces of a rectangular box region. Useful for combining with an NF2FF matrix
        // to compute the far-fields. 
        void generate_NFbox_matrix (precision x1, precision x2,
                                    precision y1, precision y2,
                                    precision z1, precision z2,
                                    int tnx, int tny, int tnz,
                                    Eigen::SparseMatrix<std::complex<precision>> &interp_mat)
        {
            Eigen::SparseMatrix<std::complex<precision>> imat_Ey_xP,imat_Ez_xP,imat_Hy_xP,imat_Hz_xP;
            Eigen::SparseMatrix<std::complex<precision>> imat_Ey_xN,imat_Ez_xN,imat_Hy_xN,imat_Hz_xN;
            Eigen::SparseMatrix<std::complex<precision>> imat_Ex_yP,imat_Ez_yP,imat_Hx_yP,imat_Hz_yP;
            Eigen::SparseMatrix<std::complex<precision>> imat_Ex_yN,imat_Ez_yN,imat_Hx_yN,imat_Hz_yN;
            Eigen::SparseMatrix<std::complex<precision>> imat_Ex_zP,imat_Ey_zP,imat_Hx_zP,imat_Hy_zP;          
            Eigen::SparseMatrix<std::complex<precision>> imat_Ex_zN,imat_Ey_zN,imat_Hx_zN,imat_Hy_zN;

            // -X
            generate_interp_matrix (x1, x1, y1, y2, z1, z2, 1, tny, tnz, COMPONENT_EY, 0, imat_Ey_xN);
            generate_interp_matrix (x1, x1, y1, y2, z1, z2, 1, tny, tnz, COMPONENT_HZ, 0, imat_Hz_xN);
            generate_interp_matrix (x1, x1, y1, y2, z1, z2, 1, tny, tnz, COMPONENT_EZ, 0, imat_Ez_xN);
            generate_interp_matrix (x1, x1, y1, y2, z1, z2, 1, tny, tnz, COMPONENT_HY, 0, imat_Hy_xN);
            // +X
            generate_interp_matrix (x2, x2, y1, y2, z1, z2, 1, tny, tnz, COMPONENT_EY, 0, imat_Ey_xP);
            generate_interp_matrix (x2, x2, y1, y2, z1, z2, 1, tny, tnz, COMPONENT_HZ, 0, imat_Hz_xP);
            generate_interp_matrix (x2, x2, y1, y2, z1, z2, 1, tny, tnz, COMPONENT_EZ, 0, imat_Ez_xP);
            generate_interp_matrix (x2, x2, y1, y2, z1, z2, 1, tny, tnz, COMPONENT_HY, 0, imat_Hy_xP);
            // -Y
            generate_interp_matrix (x1, x2, y1, y1, z1, z2, tnx, 1, tnz, COMPONENT_EZ, 0, imat_Ez_yN);
            generate_interp_matrix (x1, x2, y1, y1, z1, z2, tnx, 1, tnz, COMPONENT_HX, 0, imat_Hx_yN);
            generate_interp_matrix (x1, x2, y1, y1, z1, z2, tnx, 1, tnz, COMPONENT_EX, 0, imat_Ex_yN);
            generate_interp_matrix (x1, x2, y1, y1, z1, z2, tnx, 1, tnz, COMPONENT_HZ, 0, imat_Hz_yN); 
            // +Y
            generate_interp_matrix (x1, x2, y2, y2, z1, z2, tnx, 1, tnz, COMPONENT_EZ, 0, imat_Ez_yP);
            generate_interp_matrix (x1, x2, y2, y2, z1, z2, tnx, 1, tnz, COMPONENT_HX, 0, imat_Hx_yP);
            generate_interp_matrix (x1, x2, y2, y2, z1, z2, tnx, 1, tnz, COMPONENT_EX, 0, imat_Ex_yP);
            generate_interp_matrix (x1, x2, y2, y2, z1, z2, tnx, 1, tnz, COMPONENT_HZ, 0, imat_Hz_yP);
            // -Z
            generate_interp_matrix (x1, x2, y1, y2, z1, z1, tnx, tny, 1, COMPONENT_EX, 0, imat_Ex_zN);
            generate_interp_matrix (x1, x2, y1, y2, z1, z1, tnx, tny, 1, COMPONENT_HY, 0, imat_Hy_zN);
            generate_interp_matrix (x1, x2, y1, y2, z1, z1, tnx, tny, 1, COMPONENT_EY, 0, imat_Ey_zN);
            generate_interp_matrix (x1, x2, y1, y2, z1, z1, tnx, tny, 1, COMPONENT_HX, 0, imat_Hx_zN);
            // +Z
            generate_interp_matrix (x1, x2, y1, y2, z2, z2, tnx, tny, 1, COMPONENT_EX, 0, imat_Ex_zP);
            generate_interp_matrix (x1, x2, y1, y2, z2, z2, tnx, tny, 1, COMPONENT_HY, 0, imat_Hy_zP);
            generate_interp_matrix (x1, x2, y1, y2, z2, z2, tnx, tny, 1, COMPONENT_EY, 0, imat_Ey_zP);
            generate_interp_matrix (x1, x2, y1, y2, z2, z2, tnx, tny, 1, COMPONENT_HX, 0, imat_Hx_zP);


            // Now we need to vertically concatenate these into a single matrix.
            std::vector<Eigen::SparseMatrix<std::complex<precision>>> all_imats = {
                imat_Ey_xP,imat_Ez_xP,imat_Hy_xP,imat_Hz_xP,
                imat_Ey_xN,imat_Ez_xN,imat_Hy_xN,imat_Hz_xN,
                imat_Ex_yP,imat_Ez_yP,imat_Hx_yP,imat_Hz_yP,
                imat_Ex_yN,imat_Ez_yN,imat_Hx_yN,imat_Hz_yN,
                imat_Ex_zP,imat_Ey_zP,imat_Hx_zP,imat_Hy_zP,          
                imat_Ex_zN,imat_Ey_zN,imat_Hx_zN,imat_Hy_zN};
            
            verticalConcat_multi(all_imats, interp_mat);
        }

        void generate_interp_matrix (precision x1, precision x2, 
                                     precision y1, precision y2, 
                                     precision z1, precision z2, 
                                     int tnx, int tny, int tnz,
                                     int component, int Eonly,
                                     Eigen::SparseMatrix<std::complex<precision>> &interp_mat)
        {
            int index;
            std::vector<int> ret_x, ret_y, ret_z;
            std::vector<precision> ratio_x, ratio_y, ratio_z;
            std::vector<int> idx_x, idx_y, idx_z;
            precision step_x, step_y, step_z, cur;
            std::vector<precision> *xarr, *yarr, *zarr;
            int xidx_off, yidx_off, zidx_off;
            int xsz, ysz, zsz;
            int ehvec_off;
            precision rz1y1x1, rz1y1x2, rz1y2x1, rz1y2x2;
            precision rz2y1x1, rz2y1x2, rz2y2x1, rz2y2x2;

            if (Eonly)
                interp_mat.resize(tnx*tny*tnz, numE);
            else
                interp_mat.resize(tnx*tny*tnz, numE+numH);

            if (component == COMPONENT_EX)
            {
                xarr = &hx;
                yarr = &y;
                zarr = &z;
                xidx_off=0;
                yidx_off=-1;
                zidx_off=-1;
                ehvec_off = 0;
                xsz = Nx;
                ysz = Ny-1;
                zsz = Nz-1;
            } else if (component == COMPONENT_EY) {
                xarr = &x;
                yarr = &hy;
                zarr = &z;
                xidx_off=-1;
                yidx_off=0;
                zidx_off=-1;
                ehvec_off = numEx;
                xsz = Nx-1;
                ysz = Ny;
                zsz = Nz-1;
            } else if (component == COMPONENT_EZ) {
                xarr = &x;
                yarr = &y;
                zarr = &hz;
                xidx_off=-1;
                yidx_off=-1;
                zidx_off=0;
                ehvec_off = numEx+numEy;
                xsz = Nx-1;
                ysz = Ny-1;
                zsz = Nz;
            } else if (component == COMPONENT_HX) {
                xarr = &x;
                yarr = &hy;
                zarr = &hz;
                xidx_off=-1;
                yidx_off=0;
                zidx_off=0;              
                ehvec_off = numEx+numEy+numEz;
                xsz = Nx-1;
                ysz = Ny;
                zsz = Nz;
            } else if (component == COMPONENT_HY) {
                xarr = &hx;
                yarr = &y;
                zarr = &hz;
                xidx_off=0;
                yidx_off=-1;
                zidx_off=0;
                ehvec_off = numEx+numEy+numEz+numHx;
                xsz = Nx;
                ysz = Ny-1;
                zsz = Nz;
            } else if (component == COMPONENT_HZ) {
                xarr = &hx;
                yarr = &hy;
                zarr = &z;
                xidx_off=0;
                yidx_off=0;
                zidx_off=-1;              
                ehvec_off = numEx+numEy+numEz+numHx+numHy;
                xsz = Nx;
                ysz = Ny;
                zsz = Nz-1;
            }

            step_x = (x2-x1)/(precision)(tnx-1);
            step_y = (y2-y1)/(precision)(tny-1);
            step_z = (z2-z1)/(precision)(tnz-1);

            ret_x.resize(tnx);
            ret_y.resize(tny);
            ret_z.resize(tnz);
            idx_x.resize(tnx);
            idx_y.resize(tny);
            idx_z.resize(tnz);
            ratio_x.resize(tnx);
            ratio_y.resize(tny);
            ratio_z.resize(tnz);

            // get interpolation information for x coordinates:
            cur = x1;
            for (int ii = 0; ii < tnx; ii++)
            {
                ret_x[ii] = find_index (*xarr, cur, &ratio_x[ii], &idx_x[ii]);
                // this is because some components the 0th index is actually in the 1 cell position
                // due to the 0 cell position being a 0 boundary condition, e.g., Ey has (Nx-1) x Ny x (Nz-1)
                // elements so the first element in the y direction corresponds to y[1] rather than y[0].
                idx_x[ii] += xidx_off;
                if (!ret_x[ii])
                {
                    std::cout << "error: X position (" << ii << "): " << cur << " is out of bounds!" << std::endl;
                    exit (0);
                }
                cur += step_x;
            }
            // get interpolation information for y coordinates:
            cur = y1;
            for (int ii = 0; ii < tny; ii++)
            {
                ret_y[ii] = find_index (*yarr, cur, &ratio_y[ii], &idx_y[ii]);
                idx_y[ii] += yidx_off;

                if (!ret_y[ii])
                {
                    std::cout << "error: Y position (" << ii << "): " << cur << " is out of bounds!" << std::endl;
                    exit (0);
                }
                cur += step_y;
            }
            // get interpolation information for x coordinates:
            cur = z1;
            for (int ii = 0; ii < tnz; ii++)
            {
                ret_z[ii] = find_index (*zarr, cur, &ratio_z[ii], &idx_z[ii]);
                idx_z[ii] += zidx_off;

                if (!ret_z[ii])
                {
                    std::cout << "error: Z position (" << ii << "): " << cur << " is out of bounds!" << std::endl;
                    exit (0);
                }
                cur += step_z;
            }

            typedef Eigen::Triplet<std::complex<precision>> T;
            std::vector<T> interpMat_tL;
            interpMat_tL.reserve(8*tnx*tny*tnz);

            index=0;
            for (int kk = 0; kk < tnz; kk++)
            {
                for (int jj = 0; jj < tny; jj++)
                {
                    for (int ii = 0; ii < tnx; ii++)
                    {
                        rz1y1x1 = ratio_z[kk]*ratio_y[jj]*ratio_x[ii];
                        rz1y1x2 = ratio_z[kk]*ratio_y[jj]*(1-ratio_x[ii]);
                        rz1y2x1 = ratio_z[kk]*(1-ratio_y[jj])*ratio_x[ii];
                        rz1y2x2 = ratio_z[kk]*(1-ratio_y[jj])*(1-ratio_x[ii]);
                        rz2y1x1 = (1-ratio_z[kk])*ratio_y[jj]*ratio_x[ii];
                        rz2y1x2 = (1-ratio_z[kk])*ratio_y[jj]*(1-ratio_x[ii]);
                        rz2y2x1 = (1-ratio_z[kk])*(1-ratio_y[jj])*ratio_x[ii];
                        rz2y2x2 = (1-ratio_z[kk])*(1-ratio_y[jj])*(1-ratio_x[ii]);

                        interpMat_tL.push_back(T(index,ehvec_off+idx_x[ii]+idx_y[jj]*xsz+idx_z[kk]*xsz*ysz,(std::complex<precision>)rz1y1x1));
                        if (ret_x[ii]>1)
                            interpMat_tL.push_back(T(index,ehvec_off+(idx_x[ii]+1)+idx_y[jj]*xsz+idx_z[kk]*xsz*ysz,(std::complex<precision>)rz1y1x2));
                        if (ret_y[jj]>1)
                        {
                            interpMat_tL.push_back(T(index,ehvec_off+idx_x[ii]+(idx_y[jj]+1)*xsz+idx_z[kk]*xsz*ysz,(std::complex<precision>)rz1y2x1));
                            if (ret_x[ii]>1)
                                interpMat_tL.push_back(T(index,ehvec_off+(idx_x[ii]+1)+(idx_y[jj]+1)*xsz+idx_z[kk]*xsz*ysz,(std::complex<precision>)rz1y2x2));
                        }
                        if (ret_z[kk]>1)
                        {
                            interpMat_tL.push_back(T(index,ehvec_off+idx_x[ii]+idx_y[jj]*xsz+(idx_z[kk]+1)*xsz*ysz,(std::complex<precision>)rz2y1x1));
                            if (ret_x[ii]>1)
                                interpMat_tL.push_back(T(index,ehvec_off+(idx_x[ii]+1)+idx_y[jj]*xsz+(idx_z[kk]+1)*xsz*ysz,(std::complex<precision>)rz2y1x2));
                            if (ret_y[jj]>1)
                            {
                                interpMat_tL.push_back(T(index,ehvec_off+idx_x[ii]+(idx_y[jj]+1)*xsz+(idx_z[kk]+1)*xsz*ysz,(std::complex<precision>)rz2y2x1));
                                if (ret_x[ii]>1)
                                    interpMat_tL.push_back(T(index,ehvec_off+(idx_x[ii]+1)+(idx_y[jj]+1)*xsz+(idx_z[kk]+1)*xsz*ysz,(std::complex<precision>)rz2y2x2));
                            }    
                        }
                        index++;
                    }
                }
            }

            interp_mat.setFromTriplets(interpMat_tL.begin(), interpMat_tL.end());
        }

        // Find the inverse of the epsr matrix and store it in epsr_inv.
        void prepare_epsr_inv (void)
        {
            std::complex<precision> esx, esy, esz;
            std::complex<precision> hsx, hsy, hsz;

            // add in the PML tensors:
            for (int kk = 0; kk < Nz; kk++)
            {
                if (kk < numPMLZ)
                {
                    esz = 1.0_d + I1*sigEz[kk];
                    hsz = 1.0_d + I1*sigHz[kk];
                } else if (kk >= (Nz-numPMLZ)) {
                    if (kk==(Nz-numPMLZ))
                        esz = 1.0_d;
                    else
                        esz = 1.0_d + I1*sigEz[(Nz-1)-(kk-1)];
                    hsz = 1.0_d + I1*sigHz[(Nz-1)-kk];
                } else {
                    esz = 1.0_d;
                    hsz = 1.0_d;
                }
                for (int jj = 0; jj < Ny; jj++)
                {
                    if (jj < numPMLY)
                    {
                        esy = 1.0_d + I1*sigEy[jj];
                        hsy = 1.0_d + I1*sigHy[jj];
                    } else if (jj >= (Ny-numPMLY)) {
                        if (jj==(Ny-numPMLY))
                            esy = 1.0_d;
                        else
                            esy = 1.0_d + I1*sigEy[(Ny-1)-(jj-1)];
                        hsy = 1.0_d + I1*sigHy[(Ny-1)-jj];
                    } else {
                        esy = 1.0_d;
                        hsy = 1.0_d;
                    }    
                    for (int ii = 0; ii < Nx; ii++)
                    {
                        if (ii < numPMLX)
                        {
                            esx = 1.0_d + I1*sigEx[ii];
                            hsx = 1.0_d + I1*sigHx[ii];
                        } else if (ii >= (Nx-numPMLX)) {
                            if (ii==(Nx-numPMLX))
                                esx = 1.0_d;
                            else
                                esx = 1.0_d + I1*sigEx[(Nx-1)-(ii-1)];
                            hsx = 1.0_d + I1*sigHx[(Nx-1)-ii];
                        } else {
                            esx = 1.0_d;
                            hsx = 1.0_d;
                        }        
                        // Ex component:
                        if ((jj>0) && (kk>0))
                        {
                            epsr[0][ii+(jj-1)*Nx+(kk-1)*Nx*(Ny-1)] *= (esy*esz/hsx);
                        }
                        // Ey component:
                        if ((ii>0) && (kk>0))
                        {
                            epsr[1][(ii-1)+jj*(Nx-1)+(kk-1)*(Nx-1)*Ny] *= (esx*esz/hsy);
                        }
                        // Ez component:
                        if ((ii>0) && (jj>0))
                        {
                            epsr[2][(ii-1)+(jj-1)*(Nx-1)+kk*(Nx-1)*(Ny-1)] *= (esx*esy/hsz);
                        }
                        // Hx component:
                        if (ii>0)
                        {
                            mur[0][(ii-1)+jj*(Nx-1)+kk*(Nx-1)*Ny] *= (hsy*hsz/esx);
                        }
                        // Hy component:
                        if (jj>0)
                        {
                            mur[1][ii+(jj-1)*Nx+kk*Nx*(Ny-1)] *= (hsx*hsz/esy);
                        }
                        // Hz component:
                        if (kk>0)
                        {
                            mur[2][ii+jj*Nx+(kk-1)*Nx*Ny] *= (hsx*hsy/esz);
                        }

                    }
                }
            }

            // compute epsr^-1 from epsr vectors and
            // stuff into Eigen diagonal matrix format:
            for (int i = 0; i < numEx; i++)
            {
                if (pec[0][i])
                    epsr_inv.diagonal()[i] = 0;
                else
                    epsr_inv.diagonal()[i] = 1.0_d / epsr[0][i];
            }
            for (int i = 0; i < numEy; i++)
            {
                if (pec[1][i])
                    epsr_inv.diagonal()[i] = 0;
                else
                    epsr_inv.diagonal()[numEx+i] = 1.0_d / epsr[1][i];
            }
            for (int i = 0; i < numEz; i++)
            {
                if (pec[2][i])
                    epsr_inv.diagonal()[numEx+numEy+i] = 0;
                else 
                    epsr_inv.diagonal()[numEx+numEy+i] = 1.0_d / epsr[2][i];
            }

            for (int i = 0; i < numHx; i++)
                mur_inv.diagonal()[i] = 1.0_d / mur[0][i];
            for (int i = 0; i < numHy; i++)
                mur_inv.diagonal()[numHx+i] = 1.0_d / mur[1][i];
            for (int i = 0; i < numHz; i++)
                mur_inv.diagonal()[numHx+numHy+i] = 1.0_d / mur[2][i];
        }
     
        // Construct the DE and DH differential operator matrices
        // using triplets for the sparse matrix format
        void build_DeDh_mats_triplets (void)
        {
            // allocate space for the matrices:
            DE.resize (numHx+numHy+numHz, numEx+numEy+numEz);
            DH.resize (numEx+numEy+numEz, numHx+numHy+numHz);

            // vectors to store the triplets:
            typedef Eigen::Triplet<std::complex<precision>> T;
            std::vector<T> DE_tL, DH_tL;
            DE_tL.reserve(4*numH);
            DH_tL.reserve(4*numE);

            int index = 0;

            // Hx:
            for (int kk = 0; kk < Nz; kk++)
            {
                for (int jj = 0; jj < Ny; jj++)
                {
                    for (int ii = 0; ii < Nx-1; ii++)
                    {
                        //Dz_Ey
                        if (kk > 0)
                            DE_tL.push_back(T(index, numEx + ii+jj*(Nx-1)+(kk-1)*(Nx-1)*(Ny),-(-1.0_d / dza[kk])));
                        if (kk < (Nz-1))
                            DE_tL.push_back(T(index, numEx + ii+jj*(Nx-1)+(kk+1-1)*(Nx-1)*(Ny),-(1.0_d / dza[kk])));
                        // Dy_Ez
                        if (jj > 0)
                            DE_tL.push_back(T(index, numEx + numEy + ii+(jj-1)*(Nx-1)+kk*(Nx-1)*(Ny-1),-1.0_d / dya[jj]));
                        if (jj < (Ny-1))
                            DE_tL.push_back(T(index, numEx + numEy + ii+(jj+1-1)*(Nx-1)+kk*(Nx-1)*(Ny-1),1.0_d / dya[jj]));
                        index++;
                    }
                }
            }
            // number of Hx elements:
            std::cout << "index: " << index << "(" << numHx << ")" << std::endl;

            // Hy:            
            for (int kk = 0; kk < Nz; kk++)
            {
                for (int jj = 0; jj < Ny-1; jj++)
                {
                    for (int ii = 0; ii < Nx; ii++)
                    {
                        // Dz_Ex:
                        if (kk > 0)
                            DE_tL.push_back(T(index, ii+jj*(Nx)+(kk-1)*(Nx)*(Ny-1),-1.0_d / dza[kk]));
                        if (kk < (Nz-1))
                            DE_tL.push_back(T(index, ii+jj*(Nx)+(kk+1-1)*(Nx)*(Ny-1),1.0_d / dza[kk]));
                        // Dx_Ez:
                        if (ii > 0)
                            DE_tL.push_back(T(index, numEx + numEy + (ii-1)+jj*(Nx-1)+kk*(Nx-1)*(Ny-1),-(-1.0_d / dxa[ii])));
                        if (ii < (Nx-1))
                            DE_tL.push_back(T(index, numEx + numEy + (ii+1-1)+jj*(Nx-1)+kk*(Nx-1)*(Ny-1),-(1.0_d / dxa[ii])));
                        index++;
                    }
                }
            }
            // number of Hx+Hy elements:
            std::cout << "index: " << index << "(" << (numHx+numHy) << ")" << std::endl;

            // Hz:
            for (int kk = 0; kk < Nz-1; kk++)
            {
                for (int jj = 0; jj < Ny; jj++)
                {
                    for (int ii = 0; ii < Nx; ii++)
                    {
                        // Dy_Ex:
                        if (jj > 0)
                            DE_tL.push_back(T(index, ii+(jj-1)*(Nx)+kk*(Nx)*(Ny-1),-(-1.0_d / dya[jj])));
                        if (jj < (Ny-1))
                            DE_tL.push_back(T(index, ii+(jj+1-1)*(Nx)+kk*(Nx)*(Ny-1),-(1.0_d / dya[jj])));
                        // Dx_Ey:
                        if (ii > 0)
                            DE_tL.push_back(T(index, numEx + (ii-1)+jj*(Nx-1)+kk*(Nx-1)*(Ny),-1.0_d / dxa[ii]));
                        if (ii < (Nx-1))
                            DE_tL.push_back(T(index, numEx + (ii+1-1)+jj*(Nx-1)+kk*(Nx-1)*(Ny),1.0_d / dxa[ii]));
                        index++;
                    }
                }
            }             
            // number of Hx+Hy+Hz elements:   
            std::cout << "index: " << index << "(" << (numHx+numHy+numHz) << ")" << std::endl;

            index = 0;

            // Ex:
            for (int kk = 0; kk < Nz-1; kk++)
            {
                for (int jj = 0; jj < Ny-1; jj++)
                {
                    for (int ii = 0; ii < Nx; ii++)
                    {
                        //Dz_Hy
                        DH_tL.push_back(T(index, numHx + ii+jj*(Nx)+(kk+1-1)*(Nx)*(Ny-1),-(-1.0_d / dhza[kk])));
                        DH_tL.push_back(T(index, numHx + ii+jj*(Nx)+(kk+1)*(Nx)*(Ny-1),-(1.0_d / dhza[kk])));
                        // Dy_Hz
                        DH_tL.push_back(T(index, numHx + numHy + ii+(jj+1-1)*(Nx)+kk*(Nx)*(Ny),-1.0_d / dhya[jj]));
                        DH_tL.push_back(T(index, numHx + numHy + ii+(jj+1)*(Nx)+kk*(Nx)*(Ny),1.0_d / dhya[jj]));
                        index++;
                    }
                }
            }            
            // number of Ex elements:
            std::cout << "index: " << index << "(" << numEx << ")" << std::endl;

            // Ey:            
            for (int kk = 0; kk < Nz-1; kk++)
            {
                for (int jj = 0; jj < Ny; jj++)
                {
                    for (int ii = 0; ii < Nx-1; ii++)
                    {
                        // Dz_Hx:
                        DH_tL.push_back(T(index, ii+jj*(Nx-1)+(kk+1-1)*(Nx-1)*(Ny),-1.0_d / dhza[kk]));
                        DH_tL.push_back(T(index, ii+jj*(Nx-1)+(kk+1)*(Nx-1)*(Ny),1.0_d / dhza[kk]));
                        // Dx_Hz:
                        DH_tL.push_back(T(index, numHx + numHy + (ii+1-1)+jj*(Nx)+kk*(Nx)*(Ny),-(-1.0_d / dhxa[ii])));
                        DH_tL.push_back(T(index, numHx + numHy + (ii+1)+jj*(Nx)+kk*(Nx)*(Ny),-(1.0_d / dhxa[ii])));
                        index++;
                    }
                }
            }
            // number of Ex+Ey elements:
            std::cout << "index: " << index << "(" << (numEx+numEy) << ")" << std::endl;

            // Ez:
            for (int kk = 0; kk < Nz; kk++)
            {
                for (int jj = 0; jj < Ny-1; jj++)
                {
                    for (int ii = 0; ii < Nx-1; ii++)
                    {
                        // Dy_Hx:
                        DH_tL.push_back(T(index, ii+(jj+1-1)*(Nx-1)+kk*(Nx-1)*(Ny),-(-1.0_d / dhya[jj])));
                        DH_tL.push_back(T(index, ii+(jj+1)*(Nx-1)+kk*(Nx-1)*(Ny),-(1.0_d / dhya[jj])));
                        // Dx_Hy:
                        DH_tL.push_back(T(index, numHx + (ii+1-1)+jj*(Nx)+kk*(Nx)*(Ny-1),-1.0_d / dhxa[ii]));
                        DH_tL.push_back(T(index, numHx + (ii+1)+jj*(Nx)+kk*(Nx)*(Ny-1),1.0_d / dhxa[ii]));
                        index++;
                    }
                }
            }          
            // number of Ex+Ey+Ez elements:      
            std::cout << "index: " << index << "(" << (numEx+numEy+numEz) << ")" << std::endl;            

            // construct the sparse matrices from the triplets:
            DE.setFromTriplets(DE_tL.begin(), DE_tL.end());
            DH.setFromTriplets(DH_tL.begin(), DH_tL.end());
        }        

        // produces a matrix to compute a vector with both [E;H] fields
        // from the E-only vector, stored in EHfEmat.
        // NOTE: this currently assumes M = 0
        void setup_extractEH_matrix (void)
        {
            Eigen::SparseMatrix<std::complex<precision>> tmp1, tmp2;
            
            tmp1 = Eigen::VectorXcd::Ones(numEx+numEy+numEz).asDiagonal();
            tmp2 = mur_inv*DE/(-I1*w0*mu0);
            EHfEmat = verticalConcat(tmp1,tmp2);
        }

        // limitations: constructs PML on all 6 faces (complete open boundary)
        void setup_pml (int nPMLX, int nPMLY, int nPMLZ, precision pml_max, precision pml_pow)
        {
            precision tmp;

            numPMLX = nPMLX;
            numPMLY = nPMLY;
            numPMLZ = nPMLZ;

            sigEx.resize(numPMLX);
            sigEy.resize(numPMLY);
            sigEz.resize(numPMLZ);
            sigHx.resize(numPMLX);
            sigHy.resize(numPMLY);
            sigHz.resize(numPMLZ);
         
            for (int i = 0; i < numPMLX; i++)
            {
                tmp = (precision)(numPMLX-i)/(precision)numPMLX;
                sigEx[i] = (precision)((precision)numPMLX-((precision)i))/(precision)numPMLX;
                sigHx[i] = (precision)((precision)numPMLX-((precision)i+0.5_d))/(precision)numPMLX;
                sigEx[i]=-pml_max*pow(sigEx[i],pml_pow);
                sigHx[i]=-pml_max*pow(sigHx[i],pml_pow);
            }     
            for (int i = 0; i < numPMLY; i++)
            {
                tmp = (precision)(numPMLY-i)/(precision)numPMLY;
                sigEy[i] = (precision)((precision)numPMLY-((precision)i))/(precision)numPMLY;
                sigHy[i] = (precision)((precision)numPMLY-((precision)i+0.5_d))/(precision)numPMLY;
                sigEy[i]=-pml_max*pow(sigEy[i],pml_pow);
                sigHy[i]=-pml_max*pow(sigHy[i],pml_pow);
            }     
            for (int i = 0; i < numPMLZ; i++)
            {
                tmp = (precision)(numPMLZ-i)/(precision)numPMLZ;
                sigEz[i] = (precision)((precision)numPMLZ-((precision)i))/(precision)numPMLZ;
                sigHz[i] = (precision)((precision)numPMLZ-((precision)i+0.5_d))/(precision)numPMLZ;
                sigEz[i]=-pml_max*pow(sigEz[i],pml_pow);
                sigHz[i]=-pml_max*pow(sigHz[i],pml_pow);
            }     

        }

        // Build the system matrix A = I - epsr^-1*DH*mur^-1*DE/(w0^2*mu0*eps0)
        // after DE, DH, epsr_inv, and mur_inv have been constructed.
        void build_system_matrix (void)
        {
            Eigen::SparseMatrix<std::complex<precision>> tmp;

            tmp = (epsr_inv * DH * mur_inv * DE)/(w0*w0*mu0*eps0);
            Amat = Eigen::VectorXcd::Ones(numEx+numEy+numEz).asDiagonal();
            Amat = Amat - tmp;
        }

        ~acmefdfd ()
        {
            return;
        }
};
