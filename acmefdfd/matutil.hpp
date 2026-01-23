// Utility functions for working with Eigen sparse matrices
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

#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <string>
#include <Eigen>


// Horizontally concatenate two Eigen sparse matrices
template<typename Scalar>
Eigen::SparseMatrix<Scalar> horizontalConcat(Eigen::SparseMatrix<Scalar>& A,
                                             Eigen::SparseMatrix<Scalar>& B)
{
    // Ensure the row dimensions match.
    assert(A.rows() == B.rows());
    
    // Prepare a container for the nonzero entries.
    typedef Eigen::Triplet<Scalar> Triplet;
    std::vector<Triplet> triplets;
    triplets.reserve(A.nonZeros() + B.nonZeros());
    
    // Insert entries from A.
    for (int col = 0; col < A.outerSize(); ++col) {
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(A, col); it; ++it) {
            triplets.push_back(Triplet(it.row(), it.col(), it.value()));
        }
    }
    
    // Insert entries from B, shifting the column index by A.cols().
    for (int col = 0; col < B.outerSize(); ++col) {
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(B, col); it; ++it) {
            triplets.push_back(Triplet(it.row(), it.col() + A.cols(), it.value()));
        }
    }
    
    // Create the resulting matrix with appropriate dimensions.
    Eigen::SparseMatrix<Scalar> C(A.rows(), A.cols() + B.cols());
    C.setFromTriplets(triplets.begin(), triplets.end());
    
    return C;
}

// Vertically concatenate two Eigen sparse matrices
Eigen::SparseMatrix<std::complex<precision>> verticalConcat(Eigen::SparseMatrix<std::complex<precision>> &matrix1, Eigen::SparseMatrix<std::complex<precision>> &matrix2) 
{
    if (matrix1.cols() != matrix2.cols()) {
        throw std::invalid_argument("Matrices must have the same number of columns for vertical concatenation.");
    }

    int rows1 = matrix1.rows();
    int rows2 = matrix2.rows();
    int cols = matrix1.cols();

    Eigen::SparseMatrix<std::complex<precision>> result(rows1 + rows2, cols);

    std::vector<Eigen::Triplet<std::complex<precision>>> triplets;
    triplets.reserve(matrix1.nonZeros() + matrix2.nonZeros());

    // Copy triplets from the first matrix
    for (int k = 0; k < matrix1.outerSize(); ++k) {
        for (Eigen::SparseMatrix<std::complex<precision>>::InnerIterator it(matrix1, k); it; ++it) {
            triplets.push_back(Eigen::Triplet<std::complex<precision>>(it.row(), it.col(), it.value()));
        }
    }

    // Copy triplets from the second matrix, adjusting row indices
    for (int k = 0; k < matrix2.outerSize(); ++k) {
        for (Eigen::SparseMatrix<std::complex<precision>>::InnerIterator it(matrix2, k); it; ++it) {
            triplets.push_back(Eigen::Triplet<std::complex<precision>>(it.row() + rows1, it.col(), it.value()));
        }
    }

    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

// Vertically concatenate multiple Eigen sparse matrices provided as a vector
// of input sparse matrices
Eigen::SparseMatrix<std::complex<precision>> verticalConcat_multi(std::vector<Eigen::SparseMatrix<std::complex<precision>>> &matrix_in, Eigen::SparseMatrix<std::complex<precision>> &result) 
{
    int num_matrices, cols, tot_rows, tot_nnz;
    std::vector<int> rows(matrix_in.size());

    num_matrices = matrix_in.size();

    // check proper dimensions:
    cols = matrix_in[0].cols();
    rows[0] = matrix_in[0].rows();
    tot_rows = rows[0];
    tot_nnz = matrix_in[0].nonZeros();
    for (int i = 1; i < num_matrices; i++)
    {
        if (matrix_in[i].cols() != cols)
        {
            throw std::invalid_argument("Matrices must have the same number of columns for vertical concatenation.");
        }
        rows[i] = matrix_in[i].rows();
        tot_rows += rows[i];
        tot_nnz += matrix_in[i].nonZeros();
    }

    result.resize(tot_rows, cols);

    std::vector<Eigen::Triplet<std::complex<precision>>> triplets;
    triplets.reserve(tot_nnz);

    int row_offset;

    row_offset=0;
    for (int i = 0; i < num_matrices; i++)
    {
        // Copy triplets from the first matrix
        for (int k = 0; k < matrix_in[i].outerSize(); ++k) {
            for (Eigen::SparseMatrix<std::complex<precision>>::InnerIterator it(matrix_in[i], k); it; ++it) {
                triplets.push_back(Eigen::Triplet<std::complex<precision>>(row_offset+it.row(), it.col(), it.value()));
            }
        }
        row_offset += matrix_in[i].rows();
    }

    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

// Save a sparse matrix to a binary file with the following contents:
// 1. Number of rows (int32)
// 2. Number of columns (int32)
// 3. Number of non-zeros (int32)
// 4. Outer index array (int32 of size outerSize + 1)
// 5. Inner index array (int32 of size nnz)
// 6. Values array (complex<precision> of size nnz)
int saveSparseMatrixBinary(const Eigen::SparseMatrix<std::complex<precision>>& A, const std::string& filename)
{
    // Get dimensions and number of nonzeros.
    int rows = A.rows();
    int cols = A.cols();
    int nnz  = A.nonZeros();
    int outerSize = A.outerSize();  // For a column-major matrix, outerSize() == number of columns.

    // Open output file in binary mode.
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cout << "Error opening file " << filename << " for writing." << std::endl;
        return 0;
    }

    // Write matrix dimensions and nonzero count.
    out.write(reinterpret_cast<const char*>(&rows), sizeof(int));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(int));
    out.write(reinterpret_cast<const char*>(&nnz), sizeof(int));

    // Write the outer index array (size: outerSize + 1).
    out.write(reinterpret_cast<const char*>(A.outerIndexPtr()), sizeof(int) * (outerSize + 1));

    // Write the inner index array (size: nnz).
    out.write(reinterpret_cast<const char*>(A.innerIndexPtr()), sizeof(int) * nnz);

    // Write the values array (size: nnz). Assumes that Scalar is a POD type.
    out.write(reinterpret_cast<const char*>(A.valuePtr()), sizeof(std::complex<precision>) * nnz);

    out.close();
    return 1;
}
