// Augmented Partial Factorization (APF) solver using MUMPS
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
#include <complex>


// macros for MUMPS control arrays
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define CNTL(I) cntl[(I)-1] /* macro s.t. indices match documentation */


// wrapper for MUMPS interface
class APF_solver
{
    public:
#ifdef USE_MUMPS_SINGLE
        CMUMPS_STRUC_C id;
#else
        ZMUMPS_STRUC_C id;
#endif
        // When compiling with -DINTSIZE64, MUMPS_INT is 64-bit but MPI
        // ilp64 versions normally still require standard int for C
        // For error outputs
        int myid, ierr, error = 0;

        MUMPS_INT n;                    // number of rows in A matrix (sparse)
        MUMPS_INT8 nnz = 0;             // number of nonzeros in A 
        std::vector<MUMPS_INT> irn;     // row indices
        std::vector<MUMPS_INT> jcn;     // column indices
        std::vector<std::complex<mumps_precision>> a;         // values

        std::vector<std::complex<mumps_precision>> rhs_mumps; // APF RHS vector

        acmetime mytime;                // timekeeper

        APF_solver (void)
        {
            // MPI is not actually used, but MUMPS requires it to be initialized
            ierr = MPI_Init(NULL, NULL);
            ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

            // initialize a MUMPS instance. Use MPI_COMM_WORLD communicator
            id.comm_fortran=USE_COMM_WORLD;
            id.par=1; id.sym=0;
            id.job=JOB_INIT;
        
#ifdef USE_MUMPS_SINGLE
            cmumps_c(&id);
#else
            zmumps_c(&id);
#endif
        }

        // load the A matrix into MUMPS
        void provide_matrix (Eigen::SparseMatrix<std::complex<precision>> &A)
        {
            nnz = A.nonZeros();
            irn.resize(nnz);
            jcn.resize(nnz);
            a.resize(nnz);
    
            // grab (row, col, value) from Eigen sparse matrix
            int index=0;
            for(int i = 0; i < A.outerSize(); i++)
            {
                for(typename Eigen::SparseMatrix<std::complex<precision>>::InnerIterator it(A,i); it; ++it)
                {
                    irn[index] = it.row()+1;
                    jcn[index] = it.col()+1;
                    a[index] = (std::complex<mumps_precision>)it.value();
                    index++;
                }
            }
    
            n = A.rows();
            id.n = A.rows();
            id.nnz =nnz;
            id.irn=(MUMPS_INT*)&irn[0];
            id.jcn=(MUMPS_INT*)&jcn[0];
#ifdef USE_MUMPS_SINGLE
            id.a = (mumps_complex*)&a[0]; 
#else
            id.a = (mumps_double_complex*)&a[0]; 
#endif        
        }

        // if 1, use tree parallelism (L0-threads layer)
        // typically increases performance, but also increases memory somewhat
        void use_L0_threads (int use)
        {
            id.ICNTL(48) = use&1;
        }

        // set number of OpenMP threads to use for factorization
        void set_num_threads (int num_threads)
        {
            id.ICNTL(16) = num_threads;
        }

        // solve for the Schur complement with dimensions size_schur^2 
        // and store in S
        void solve_apf (Eigen::MatrixXcd &S, int size_schur)
        {
            std::vector<int> schur_inds(size_schur);

            std::vector<std::complex<mumps_precision>> schur_scratch(size_schur*size_schur);

            id.size_schur = size_schur;
            for (int i = 0; i < size_schur; i++)
                schur_inds[i] = n-size_schur+i+1;
            id.listvar_schur = (MUMPS_INT*)&schur_inds[0];
            id.nprow = 1;
            id.npcol = 1;
            id.mblock = 100;
            id.nblock = 100;

#ifdef USE_MUMPS_SINGLE
            id.schur = (mumps_complex*)&schur_scratch[0]; 
#else
            id.schur = (mumps_double_complex*)&schur_scratch[0]; 
#endif
            id.schur_lld = id.size_schur;
            id.ICNTL(19) = 3; // centralized Schur factors by columns
            id.ICNTL(31) = 1; // discard factors
            // error and global info to stdout, no diagnostic prints
            id.ICNTL(1)=6; id.ICNTL(2)=0; id.ICNTL(3)=6; id.ICNTL(4)=2;
            id.ICNTL(7) = 5; // use Metis to compute ordering
            id.ICNTL(14) = 35; // max % increase in working memory estimate

            // do analysis:
            id.job = 1;
            mytime.tik();
#ifdef USE_MUMPS_SINGLE
            cmumps_c(&id);
#else
            zmumps_c(&id);
#endif

            // if error, print INFOG error information
            if (id.infog[0]<0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                myid, id.infog[0], id.infog[1]);
                error = 1;
                return;
            }
            // print elapsed time
            double diff = mytime.tok();
            std::cout << "MUMPS APF analyze time elapsed: " << diff << std::endl;  

            // do factorization:
            id.job = 2;
            mytime.tik();
#ifdef USE_MUMPS_SINGLE
            cmumps_c(&id);
#else
            zmumps_c(&id);
#endif
            // if error, print INFOG error information
            if (id.infog[0]<0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                myid, id.infog[0], id.infog[1]);
                error = 1;
                return;
            }
            // print elapsed time
            diff = mytime.tok();
            std::cout << "MUMPS APF factorize time elapsed: " << diff << std::endl;

            // write Schur data out to the provided matrix:
            S.resize(size_schur,size_schur);
            int index=0;
            for (int jj = 0; jj < size_schur; jj++)
            {
                for (int ii = 0; ii < size_schur; ii++)
                {
                    S(ii,jj) = -(std::complex<precision>)schur_scratch[index]; //mumps_schur[index];
                    index++;
                }
            }
        }

        // on return, res vector holds result.
        void solve_rhs (Eigen::VectorXcd &rhs, Eigen::VectorXcd &res)
        {
            rhs_mumps.resize(rhs.rows());
            res.resize(rhs.rows());

            for (int i = 0; i < rhs.rows(); i++)
                rhs_mumps[i] = (std::complex<mumps_precision>)rhs(i);

#ifdef USE_MUMPS_SINGLE
            id.rhs = (mumps_complex*)&rhs_mumps[0];
#else 
            id.rhs = (mumps_double_complex*)&rhs_mumps[0];
#endif        

            // error and global info to stdout, no diagnostic prints
            id.ICNTL(1)=6; id.ICNTL(2)=0; id.ICNTL(3)=6; id.ICNTL(4)=2;
            id.ICNTL(7) = 5; // use Metis to compute ordering
            id.ICNTL(14) = 35; // max % increase in working memory estimate 

            error=0;
            mytime.tik();
            // Call the MUMPS package (analyze, factorization, and solve)
            id.job=6;
#ifdef USE_MUMPS_SINGLE
                cmumps_c(&id);
#else
                zmumps_c(&id);
#endif
            // if error, print INFOG error information
            if (id.infog[0]<0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                myid, id.infog[0], id.infog[1]);
                error = 1;
            }
            // print elapsed time
            double diff = mytime.tok();
            std::cout << "MUMPS solve time elapsed: " << diff << std::endl; 
            for (int i = 0; i < rhs.rows(); i++)
                res(i) = (std::complex<precision>)rhs_mumps[i];          
        }

        ~APF_solver ()
        {
            // terminate the MUMPS instance
            id.job=JOB_END;
#ifdef USE_MUMPS_SINGLE
            cmumps_c(&id);
 #else
            zmumps_c(&id);
#endif                   
        }
};
