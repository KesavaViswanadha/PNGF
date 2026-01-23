// Finite Difference Based Modesolver
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

class acmemode
{
    public:
        // inputs
        int Nx, Ny;                 // number of Ez cells in x and y
        double f0, w0, lambda0;     // frequency
        double dx, dy;              // cell size

        int num_i_Hx, num_i_Hy;     // number of components on Yee cell grid
        int num_i_Ez, num_i_Hz;
        int num_j_Hx, num_j_Hy;
        int num_j_Ez, num_j_Hz;

        std::complex<precision> betasq; // propagation constants (to be solved)

        // Build the whole L block matrix: L = [L11 L12; L21 L22].
        Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> L;

        // Nx and Ny correspond to the number of Ez cells.
        // Outside boundary is terminated in PEC BCs.
        acmemode (int mNx, int mNy, precision mf0, precision mdx, precision mdy)
        {
            Nx = mNx-1;
            Ny = mNy-1;

            f0 = mf0;
            w0 = 2*M_PI*f0;

            lambda0 = c0/f0;
            dx = lambda0/50.0;
            dy = lambda0/50.0;

            // Enumerate the number of components on Yee cell grid:
            // . = Ez, x = Hz, - = Hx, | = Hy
            // (Cells on edge are not unknowns since they are forced to 0 to satisfy PEC BCs.)
            //
            // . | . | .
            // - x - x -
            // . | . | . 
            // - x - x -
            // . | . | . 
            //
            num_i_Hx = Nx-1;
            num_i_Hy = Nx;
            num_j_Hx = Ny;
            num_j_Hy = Ny-1;
            num_i_Ez = Nx-1;
            num_j_Ez = Ny-1;
            num_i_Hz = Nx;
            num_j_Hz = Ny;    
        }

        // construct the finite difference matrices needed to find the 
        // propagation constants of the modes
        int build_operator (void)
        {
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> Dx_Hy;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> Dy_Hx;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> Dy_Ez;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> Dx_Ez;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> Dx_Hx;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> Dy_Hy;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> Dy_Hz;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> Dx_Hz;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> epsr_ez;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> epsr_hx;
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> epsr_hy;

            // Allocate the differentiation matrices:
            Dx_Hy.resize(num_i_Ez*num_j_Ez, num_i_Hy*num_j_Hy);
            Dy_Hx.resize(num_i_Ez*num_j_Ez, num_i_Hx*num_j_Hx);

            Dy_Ez.resize(num_i_Hx*num_j_Hx, num_i_Ez*num_j_Ez);
            Dx_Ez.resize(num_i_Hy*num_j_Hy, num_i_Ez*num_j_Ez);

            Dx_Hx.resize(num_i_Hz*num_j_Hz, num_i_Hx*num_j_Hx);
            Dy_Hy.resize(num_i_Hz*num_j_Hz, num_i_Hy*num_j_Hy);

            Dy_Hz.resize(num_i_Hy*num_j_Hy, num_i_Hz*num_j_Hz);
            Dx_Hz.resize(num_i_Hx*num_j_Hx, num_i_Hz*num_j_Hz);

            // Allocate the epsr matrices:
            epsr_ez.resize(num_i_Ez*num_j_Ez, num_i_Ez*num_j_Ez);
            epsr_hx.resize(num_i_Hx*num_j_Hx, num_i_Hx*num_j_Hx);
            epsr_hy.resize(num_i_Hy*num_j_Hy, num_i_Hy*num_j_Hy);

            // compute the finite differences and insert into matrices
            int idx=0;
            for (int jj = 0; jj < num_j_Ez; jj++)
            {
                for (int ii = 0; ii < num_i_Ez; ii++)
                {
                    Dx_Hy.insert(idx, jj*num_i_Hy + ii) = -1.0/dx;
                    Dx_Hy.insert(idx, jj*num_i_Hy + ii + 1) = 1.0/dx;

                    Dy_Hx.insert(idx, jj*num_i_Hx + ii) = -1.0/dy;
                    Dy_Hx.insert(idx, (jj+1)*num_i_Hx + ii) = 1.0/dy;

                    idx++;
                }
            }

            idx=0;
            for (int jj = 0; jj < num_j_Hz; jj++)
            {
                for (int ii = 0; ii < num_i_Hz; ii++)
                {
                    if (ii<(num_i_Hz-1))
                        Dx_Hx.insert(idx, jj*num_i_Hx + ii) = 1.0/dx;
                    if (ii>0)
                        Dx_Hx.insert(idx, jj*num_i_Hx + ii - 1) = -1.0/dx;

                    if (jj<(num_j_Hz-1))
                        Dy_Hy.insert(idx, jj*num_i_Hy + ii) = 1.0/dy;
                    if (jj>0)
                        Dy_Hy.insert(idx, (jj-1)*num_i_Hy + ii) = -1.0/dy;

                    idx++;
                }
            }    

            idx=0;
            for (int jj = 0; jj < num_j_Hx; jj++)
            {
                for (int ii = 0; ii < num_i_Hx; ii++)
                {
                    if (jj<(num_j_Hx-1))
                        Dy_Ez.insert(idx, jj*num_i_Ez + ii) = 1.0/dy;
                    if (jj>0)
                        Dy_Ez.insert(idx, (jj-1)*num_i_Ez + ii) = -1.0/dy;

                    idx++;
                }
            }    

            idx=0;
            for (int jj = 0; jj < num_j_Hy; jj++)
            {
                for (int ii = 0; ii < num_i_Hy; ii++)
                {
                    if (ii<(num_i_Hy-1))
                        Dx_Ez.insert(idx, jj*num_i_Ez + ii) = 1.0/dx;
                    if (ii>0)
                        Dx_Ez.insert(idx, jj*num_i_Ez + ii - 1) = -1.0/dx;

                    idx++;
                }
            }    

            idx=0;
            for (int jj = 0; jj < num_j_Hx; jj++)
            {
                for (int ii = 0; ii < num_i_Hx; ii++)
                {
                    Dx_Hz.insert(idx, jj*num_i_Hz + ii + 1) = 1.0/dx;
                    Dx_Hz.insert(idx, jj*num_i_Hz + ii) = -1.0/dx;
                    idx++;
                }
            }    

            idx=0;
            for (int jj = 0; jj < num_j_Hy; jj++)
            {
                for (int ii = 0; ii < num_i_Hy; ii++)
                {
                    Dy_Hz.insert(idx, (jj+1)*num_i_Hz + ii) = 1.0/dy;
                    Dy_Hz.insert(idx, jj*num_i_Hz + ii) = -1.0/dy;
                    idx++;
                }
            }    

            // epsilon matrices are diagonal:
            epsr_ez = (1.0/eps0)*Eigen::VectorXcd::Ones(num_i_Ez*num_j_Ez).asDiagonal();
            epsr_hx = eps0*Eigen::VectorXcd::Ones(num_i_Hx*num_j_Hx).asDiagonal();
            epsr_hy = eps0*Eigen::VectorXcd::Ones(num_i_Hy*num_j_Hy).asDiagonal();

            // L is a linear operator such that
            //      L [hx; hy] = beta^2 [hx; hy]
            // and L = [L11 L12; L21 L22], where L11 - L22 are given by the
            // waveguide mode equations
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> L11, L12, L21, L22, tmp, tmp1, tmp2;

            tmp1 = Dy_Ez*epsr_ez;
            tmp2 = -epsr_hx*tmp1;

            tmp1 = tmp2*(-Dy_Hx);
            tmp = Dx_Hz*Dx_Hx;
            L11 = (tmp1 + tmp + w0*w0*epsr_hx*mu0).pruned();

            tmp1 = tmp2*(Dx_Hy);
            tmp = Dx_Hz*Dy_Hy;
            L12 = (tmp1 + tmp).pruned();

            tmp1 = Dx_Ez*epsr_ez;
            tmp2 = epsr_hy*tmp1;

            tmp1 = tmp2*(-Dy_Hx);
            tmp = Dy_Hz*Dx_Hx;
            L21 = (tmp1 + tmp).pruned();

            tmp1 = tmp2*(Dx_Hy);
            tmp = Dy_Hz*Dy_Hy;
            L22 = (tmp1 + tmp + w0*w0*epsr_hy*mu0).pruned();

            // concatenate the four matrices together into L:
            L.resize(L11.rows() + L22.rows(), L11.cols() + L12.cols());
            L.resizeNonZeros (L11.nonZeros()+L12.nonZeros()+L21.nonZeros()+L22.nonZeros());

            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>::StorageIndex i = 0;

            // horizontally combine L11 and L12:
            for (Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>::StorageIndex row = 0; row < L11.rows(); row++)
            {
                L.outerIndexPtr()[row] = i;

                for (Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>::StorageIndex j = L11.outerIndexPtr()[row]; j < L11.outerIndexPtr()[row + 1]; j++) //, i++)
                {
                    L.innerIndexPtr()[i] = L11.innerIndexPtr()[j];
                    L.valuePtr()[i] = L11.valuePtr()[j];
                    i++;
                }
                for (Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>::StorageIndex j = L12.outerIndexPtr()[row]; j < L12.outerIndexPtr()[row + 1]; j++) //, i++)
                {
                    L.innerIndexPtr()[i] = L11.cols() + L12.innerIndexPtr()[j];
                    L.valuePtr()[i] = L12.valuePtr()[j];
                    i++;
                }
            }
            // now combine and stack L21 and L22:
            for (Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>::StorageIndex row = 0; row < L22.rows(); row++)
            {
                L.outerIndexPtr()[L11.rows()+row] = i;

                for (Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>::StorageIndex j = L21.outerIndexPtr()[row]; j < L21.outerIndexPtr()[row + 1]; j++) //, i++)
                {
                    L.innerIndexPtr()[i] = L21.innerIndexPtr()[j];
                    L.valuePtr()[i] = L21.valuePtr()[j];
                    i++;
                }
                for (Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>::StorageIndex j = L22.outerIndexPtr()[row]; j < L22.outerIndexPtr()[row + 1]; j++) //, i++)
                {
                    L.innerIndexPtr()[i] = L21.cols() + L22.innerIndexPtr()[j];
                    L.valuePtr()[i] = L22.valuePtr()[j];
                    i++;
                }
            }
            L.outerIndexPtr()[L11.rows()+L22.rows()] = i;
            assert(L.isCompressed());   
            return 1;     
        }

        // run the mode solver once the operator is built
        // beta_guess is an initial guess for the propagation constant
        // the eigenvalue solver solves the linear system for beta and the
        // transverse magnetic field components hx and hy
        std::complex<precision> solve (std::complex<precision> beta_guess, precision reltol = 1.0E-10, int maxiter = 32)
        {
            Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> tmp;

            tmp = (-beta_guess)*Eigen::VectorXcd::Ones(L.rows()).asDiagonal();
            tmp = tmp + L;

            // eigenvalue solver:
            Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>, Eigen::COLAMDOrdering<int> >   solver;
            // fill A and b;
            // Compute the ordering permutation vector from the structural pattern of A
            solver.analyzePattern(tmp); 
            // Compute the numerical factorization 
            solver.factorize(tmp); 

            Eigen::VectorXcd b = Eigen::VectorXcd::Ones(L.cols()), x, b_prev, err_vec;

            //Use the factors to solve the linear system 

            b_prev = b;
            int iter;
            precision err;
            for (iter = 0; iter < maxiter; iter++)
            {
                x = solver.solve(b);
                b = x/x.norm();
                err_vec = b-b_prev;
                err = err_vec.norm() / b.norm();
                b_prev = b;
                if (err < reltol)
                {
                    std::cout << "converged with tol: " << err << std::endl; 
                    break;
                }
                //cout << "err:\t" << err << endl;
            }
            if (iter == maxiter)
            {
                std::cout << "maxiter reached at tol: " << err << std::endl;
            }

            std::complex<double> tcd;

            betasq = (b.dot(L * b));
            tcd = b.squaredNorm();
            betasq /= tcd;
            return std::sqrt(betasq);
        }
};
