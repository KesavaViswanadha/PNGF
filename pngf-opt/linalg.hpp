// Linear Algebra helper class. Wraps BLAS/Lapack functions:
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

class LinAlg
{
public:
    static void matrixVectorMult (Eigen::Matrix<std::complex<precision>,Dynamic,Dynamic> &A, Eigen::Vector<std::complex<precision>,Dynamic> &v,  Eigen::Vector<std::complex<precision>,Dynamic> &res)
    {
        std::complex<precision> alpha = 1.0;
        std::complex<precision> beta  = 0.0;   
#ifdef USE_SINGLE_PRECISION     
        cblas_cgemv(CblasColMajor, CblasNoTrans, A.rows(), A.cols(), &alpha,
            A.data(), A.rows(), v.data(), 1, &beta, res.data(), 1);
#else
        cblas_zgemv(CblasColMajor, CblasNoTrans, A.rows(), A.cols(), &alpha,
        A.data(), A.rows(), v.data(), 1, &beta, res.data(), 1);
#endif
    }

    static int invertMatrix (Eigen::Matrix<std::complex<precision>,Dynamic,Dynamic> &A)
    {
        int info, m, n;
        std::complex<precision> alpha = 1.0;
        std::complex<precision> beta  = 0.0;   
        std::complex<precision> work_query;
        std::vector<std::complex<precision>> work;

        int lwork = -1;

        m = A.rows();
        n = A.cols();
        std::cout << "size: " << m << ", " << n << std::endl;
        std::vector<int> ipiv(n);

        // Perform LU factorization using zgetrf_
#ifdef USE_SINGLE_PRECISION  
        cgetrf_(&n, &n, reinterpret_cast<lp_complex*>(A.data()), &n, &ipiv[0], &info);
#else
        zgetrf_(&n, &n, reinterpret_cast<lp_complex*>(A.data()), &n, &ipiv[0], &info);
#endif
        if (info != 0) {
            std::cout << "Error in zgetrf_: " << info << std::endl;
            return 0;
        }

#ifdef USE_SINGLE_PRECISION     
        cgetri_ (&n, reinterpret_cast<lp_complex*>(A.data()), &n, &ipiv[0], reinterpret_cast<lp_complex*>(&work_query), &lwork, &info);
#else
        zgetri_ (&n, reinterpret_cast<lp_complex*>(A.data()), &n, &ipiv[0], reinterpret_cast<lp_complex*>(&work_query), &lwork, &info);
#endif
        if (info != 0) {
            std::cout << "Error in zgetri_ workspace query: " << info << std::endl;
            return 0;
        }
        // The optimal workspace size is returned in work_query.r (real part).
        lwork = work_query.real();
        work.resize(lwork);

        // Compute the inverse of the matrix using zgetri_
#ifdef USE_SINGLE_PRECISION     
        cgetri_(&n, reinterpret_cast<lp_complex*>(A.data()), &n, &ipiv[0], reinterpret_cast<lp_complex*>(&work[0]), &lwork, &info);
#else  
        zgetri_(&n, reinterpret_cast<lp_complex*>(A.data()), &n, &ipiv[0], reinterpret_cast<lp_complex*>(&work[0]), &lwork, &info);
#endif
        std::cout << "info: " << info << std::endl;
        if (info != 0) {
            std::cout << "Error in zgetri_: " << info << std::endl;
            return 0;
        }
        return 1;
    }    
};