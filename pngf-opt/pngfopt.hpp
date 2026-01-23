// This class implements the PNGF method for inverse design 
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

class PNGF_Optimizer
{
public:
    Eigen::Matrix<std::complex<precision>,Dynamic,Dynamic> G; // the PNGF matrix
    Eigen::Matrix<std::complex<precision>,Dynamic,Dynamic> C; // the C matrix
    Eigen::Matrix<std::complex<precision>,Dynamic,Dynamic> Gobj; // the objective function matrix
    Eigen::Matrix<std::complex<precision>,Dynamic,Dynamic> R; // the R matrix
    Eigen::Matrix<std::complex<precision>,Dynamic,1> Einc; // the objective function matrix
    Eigen::Matrix<std::complex<precision>,Dynamic,1> PEinc, S, Snew; // the objective function matrix    
    std::vector<complex<precision>> VinvCU, RU, CinvU, VCinv, V, VSnew; // the objective function matrix    
    std::vector<int> tile_map;
    std::vector<int> yee_grid, diff_P, P_indices; // this is effectively the "P" matrix.
    std::vector<std::complex<precision>> work;

    int map_xsz, map_ysz; // number of tiles in x/y
    int map_sz; // total number of tiles.
    // number of Yee cells comprising a tile on the map:
    int tiles_yee_x, tiles_yee_y;
    // size of corresponding yee grid:
    int yee_x_xsz, yee_x_ysz, yee_x_sz;
    int yee_y_xsz, yee_y_ysz, yee_y_sz;
    int yee_sz;
    std::vector<int> port_indices;
    int port_width, port_height;
    std::complex<precision> Vinc, LN1inc, LN2inc;
    Eigen::Matrix<std::complex<precision>,Dynamic,1> Eobj, newEobj;
    std::vector<int> ipiv; // pivot scratch space
    int num_yee_flipped;
    precision dx, dy, dz;
    int last_flipped[10];
    int num_flipped_tiles;
    int num_eobj;
    int max_flip_sz;

    PNGF_Optimizer()
    {
    }

    PNGF_Optimizer (const char *filename, int num_tiles_x, int num_tiles_y,
        int num_cells_per_tile_x, int num_cells_per_tile_y,
        int tnum_eobj, int max_flips,
        precision mdx, precision mdy, precision mdz)
    {
        init(filename, num_tiles_x, num_tiles_y,
        num_cells_per_tile_x, num_cells_per_tile_y,
        tnum_eobj, max_flips, mdx, mdy, mdz);
    }    

    void init (const char *filename, int num_tiles_x, int num_tiles_y,
                    int num_cells_per_tile_x, int num_cells_per_tile_y,
                    int tnum_eobj, int max_flips,
                    precision mdx, precision mdy, precision mdz)
    {
        dx = mdx;
        dy = mdy;
        dz = mdz;

        num_eobj = tnum_eobj;

        map_xsz = num_tiles_x;
        map_ysz = num_tiles_y;
        map_sz = map_xsz*map_ysz;

        max_flip_sz = (3*4)*2*max_flips;

        tiles_yee_x = num_cells_per_tile_x;
        tiles_yee_y = num_cells_per_tile_y;

        yee_x_xsz = map_xsz*tiles_yee_x;
        yee_x_ysz = map_ysz*tiles_yee_y+1;
        yee_y_xsz = map_xsz*tiles_yee_x+1;
        yee_y_ysz = map_ysz*tiles_yee_y;
        
        yee_x_sz = yee_x_xsz*yee_x_ysz;
        yee_y_sz = yee_y_xsz*yee_y_ysz;

        yee_sz = yee_x_sz + yee_y_sz;

        tile_map.resize(num_tiles_x*num_tiles_y);
        yee_grid.resize(yee_sz);
        diff_P.resize(yee_sz);
        P_indices.resize(yee_sz);

        G.resize(yee_sz,yee_sz);
        C.resize(yee_sz,yee_sz);
        VinvCU.resize(max_flip_sz*max_flip_sz); // max possible size.
        CinvU.resize(max_flip_sz*yee_sz);
        VCinv.resize(max_flip_sz*yee_sz);
        V.resize(max_flip_sz*yee_sz);
        RU.resize(num_eobj*yee_sz);
        ipiv.resize(yee_sz);
        S.resize(yee_sz);
        Snew.resize(yee_sz);
        VSnew.resize(yee_sz);
        Gobj.resize(num_eobj,yee_sz);
        Eobj.resize(num_eobj,1);
        newEobj.resize(num_eobj,1);
        R.resize(num_eobj,yee_sz);
        Einc.resize(yee_sz,1);
        PEinc.resize(yee_sz,1);

        work.resize(100*yee_sz);

        // read in the G matrix:
        FILE *ifil = fopen(filename, "rb");
        if (!ifil)
        {
            std::cout << "Error file " << filename << " not found!" << std::endl;
            return;
        }
#ifdef USE_SINGLE_PRECISION
#if USE_ORIG_GMAT
        std::vector<std::complex<double>> tmp_buf(yee_sz*yee_sz);
        fread (&tmp_buf[0], 1, sizeof(std::complex<double>)*yee_sz*yee_sz, ifil);
        #pragma omp parallel for
        for (int ii = 0; ii < yee_sz*yee_sz; ii++)
            G.data()[ii] = (std::complex<precision>)tmp_buf[ii];
#else
            // The voltage on the source port is always assumed to be a desired observable so it is counted in
            // num_eobj, but we already have it from the main G matrix.
            std::vector<std::complex<double>> tmp_buf((yee_sz+num_eobj-1)*(yee_sz+num_eobj-1));
            fread (&tmp_buf[0], 1, sizeof(std::complex<double>)*(yee_sz+num_eobj-1)*(yee_sz+num_eobj-1), ifil);
            #pragma omp parallel for
            for (int ii = 0; ii < yee_sz*yee_sz; ii++)
            {
                int col = ii / yee_sz;
                int row = ii - col * yee_sz;
    
                G.data()[ii] = (std::complex<precision>)tmp_buf[row+(yee_sz+num_eobj-1)*col];
            }
            #pragma omp parallel for
            for (int ii = 0; ii < yee_sz*(num_eobj-1); ii++)
            {
                int col = ii / (num_eobj-1);
                int row = ii - col * (num_eobj-1);
    
                Gobj(1+row,col) = (std::complex<precision>)tmp_buf[(yee_sz+row)+(yee_sz+num_eobj-1)*col];
            }
#endif            
#else
#if USE_ORIG_GMAT
        fread (G.data(), 1, sizeof(std::complex<precision>)*yee_sz*yee_sz, ifil);
#else
        // The voltage on the source port is always assumed to be a desired observable so it is counted in
        // num_eobj, but we already have it from the main G matrix.
        std::vector<std::complex<double>> tmp_buf((yee_sz+num_eobj-1)*(yee_sz+num_eobj-1));
        fread (&tmp_buf[0], 1, sizeof(std::complex<double>)*(yee_sz+num_eobj-1)*(yee_sz+num_eobj-1), ifil);
        #pragma omp parallel for
        for (int ii = 0; ii < yee_sz*yee_sz; ii++)
        {
            int col = ii / yee_sz;
            int row = ii - col * yee_sz;

            G.data()[ii] = tmp_buf[row+(yee_sz+num_eobj-1)*col];
        }
        #pragma omp parallel for
        for (int ii = 0; ii < yee_sz*(num_eobj-1); ii++)
        {
            int col = ii / (num_eobj-1);
            int row = ii - col * (num_eobj-1);

            Gobj(1+row,col) = tmp_buf[(yee_sz+row)+(yee_sz+num_eobj-1)*col];
        }
#endif

#endif
        fclose (ifil);

        // defaults:
        port_width = tiles_yee_x;
        port_height = 2;

        return;
    }

    // Compute the initial C matrix from G based
    // on the tilemap (P matrix) and G:
    // C = (I - P) + PG from the paper.
    void compute_C_matrix (void)
    {
        int base;
#ifdef ROW_MAJOR_FORMAT
        #pragma omp parallel for
        for (int jj = 0; jj < yee_sz; jj++)
        {
            base = jj*yee_sz;
            if (yee_grid[jj]==1)
            {
                for (int ii = 0; ii < yee_sz; ii++)
                    C.data()[base+ii] = G.data()[base+ii];
            } else {
                for (int ii = 0; ii < yee_sz; ii++)
                    C.data()[base+ii] = 0.0_d;
                C.data()[base+jj] = 1.0_d;
            }    
        }
#else
        #pragma omp parallel for
        for (int ii = 0; ii < yee_sz; ii++)
        {
            base = ii*yee_sz;
            for (int jj = 0; jj < yee_sz; jj++)
            {
                if (yee_grid[jj]==1)
                {
                    C.data()[base+jj] = G.data()[base+jj];
                } else {
                    if (ii == jj)
                        C.data()[base+jj] = 1.0_d;
                    else
                        C.data()[base+jj] = 0.0_d;    
                }
            }
        }
#endif        
    }

    // Compute R = Gobj*Cinv matrix from the paper.
    // if compute_Eobj = 1, this also computes Eobj by doing R*PEinc.
    void compute_R_matrix (int compute_Eobj)
    {
        // should be called when "C" corresponds to inv(C).
        R = Gobj*C;
        if (compute_Eobj)
        {
            S = C*PEinc;
            Eobj = R*PEinc;
            newEobj = Eobj;
        } else {
            S = C*PEinc;
            Eobj = newEobj;
        }
    }

    // Use Woodbury matrix formula to update Cinv matrix based on the last tile(s) flipped
    // by the compute_updated_Eobj function. This function should NOT be called unless if
    // compute_updated_Eobj has been called at least once since the last time Cinv was
    // updated, as the dP matrix will not be valid otherwise:
    void update_C_inv (void)
    {
        int lwork;
 
        // compute Cinv*U:
        #pragma omp parallel for        
        for (int ii = 0; ii < num_yee_flipped; ii++)
        {
            int base = ii*yee_sz;
            int C_base = P_indices[ii]*yee_sz;
            if (diff_P[P_indices[ii]]>0)
            {
                for (int jj = 0; jj < yee_sz; jj++)
                    CinvU.data()[base+jj] = C.data()[C_base+jj];
            } else {
                for (int jj = 0; jj < yee_sz; jj++)
                    CinvU.data()[base+jj] = -C.data()[C_base+jj];
            }
        }

        // now compute inv(I+VinvCU)*V:
        int nrhs = 1; // number of right-hand sides (we have a single vector b)
        int lda = num_yee_flipped;  // leading dimension of A
        int ldb = num_yee_flipped;  // leading dimension of b (or x)
        int info = 0; // output status

        // the matrix has already been presumably factorized from the last call to
        // compute_updated_Eobj, so compute its inverse with ?getri:
        lwork = 100*num_yee_flipped;
#ifdef USE_SINGLE_PRECISION
        cgetri_(&num_yee_flipped, reinterpret_cast<lp_complex*>(&VinvCU[0]), &num_yee_flipped, &ipiv[0],
                reinterpret_cast<lp_complex*>(&work[0]), &lwork, &info);
#else
        zgetri_(&num_yee_flipped, reinterpret_cast<lp_complex*>(&VinvCU[0]), &num_yee_flipped, &ipiv[0], 
                reinterpret_cast<lp_complex*>(&work[0]), &lwork, &info);
#endif

        // compute V*Cinv:
        #pragma omp parallel for        
        for (int ii = 0; ii < num_yee_flipped*yee_sz; ii++)
        {
            int col = ii/num_yee_flipped;
            int row = ii - col*num_yee_flipped;

            int base_V = row*yee_sz;
            int base_C = col*yee_sz;
            std::complex<precision> ctmp;
            ctmp = 0;
            for (int i = 0; i < yee_sz; i++)
                ctmp += V.data()[base_V+i] * C.data()[base_C+i];
            VCinv.data()[ii] = ctmp;                
        }

        std::complex<precision> alpha = {1.0, 0.0};
        std::complex<precision> beta  = {0, 0.0};

        // now do: VinVCU * VCinv:
#ifdef USE_SINGLE_PRECISION
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            num_yee_flipped,    // m: number of rows of A and C.
            yee_sz,    // n: number of columns of B and C.
            num_yee_flipped,    // k: number of columns of A and rows of B.
            &alpha,
            reinterpret_cast<lp_complex*>(VinvCU.data()), num_yee_flipped,
            reinterpret_cast<lp_complex*>(VCinv.data()), num_yee_flipped,
            &beta,
            reinterpret_cast<lp_complex*>(V.data()), num_yee_flipped); // clobber V
#else
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            num_yee_flipped,    // m: number of rows of A and C.
            yee_sz,    // n: number of columns of B and C.
            num_yee_flipped,    // k: number of columns of A and rows of B.
            &alpha,
            reinterpret_cast<lp_complex*>(VinvCU.data()), num_yee_flipped,
            reinterpret_cast<lp_complex*>(VCinv.data()), num_yee_flipped,
            &beta,
            reinterpret_cast<lp_complex*>(V.data()), num_yee_flipped); // clobber V
#endif
        alpha = {-1.0, 0.0};
        beta  = {1.0, 0.0};

#ifdef USE_SINGLE_PRECISION
        // now do: C = C - Cinv*U * (VinvCU * VCinv)
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            yee_sz,   // m: number of rows of A and C.
            yee_sz,   // n: number of columns of B and C.
            num_yee_flipped,   // k: number of columns of A and rows of B.
            &alpha,
            reinterpret_cast<lp_complex*>(CinvU.data()), yee_sz,
            reinterpret_cast<lp_complex*>(V.data()), num_yee_flipped,
            &beta,
            reinterpret_cast<lp_complex*>(C.data()), yee_sz);
#else
        // now do: C = C - Cinv*U * (VinvCU * VCinv)
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            yee_sz,   // m: number of rows of A and C.
            yee_sz,   // n: number of columns of B and C.
            num_yee_flipped,   // k: number of columns of A and rows of B.
            &alpha,
            reinterpret_cast<lp_complex*>(CinvU.data()), yee_sz,
            reinterpret_cast<lp_complex*>(V.data()), num_yee_flipped,
            &beta,
            reinterpret_cast<lp_complex*>(C.data()), yee_sz);
#endif
        // Update the tile map to reflect the permanent flipped tile changes:
        for (int i = 0; i < num_flipped_tiles; i++)
            tile_map[last_flipped[i]] ^= 1;
        // P_n = P_{n-1} + dP
        #pragma omp parallel for
        for (int ii = 0; ii < yee_sz; ii++)
            yee_grid[ii] += diff_P[ii];
    }

    // Compute change in Eobj due to flipping one or more tiles.
    // Note: This is an O(Nopt) operation as it does not need to compute
    // the full updated matrix Cinv. If the change is favorable, call the
    // update_C_inv function to update Cinv.
    // Inputs: tile_idx: Array of integers containing the indices of the tiles to be flipped
    //         num_tiles_to_flip: Size of tile_idx array (should be minimum 1)
    // Output: The first entry in Eobj is returned. 
    // Note: Eobj can have more than one entry for multi-objective function optimizations
    // (e.g., optimizign S11 and gain of an antenna), but to read the other entries, the
    // newEobj object of this class should be accessed directly.
    std::complex<precision> compute_updated_Eobj (int *tile_idx, int num_tiles_to_flip)
    {
        // temporarily flip the tile:
        num_flipped_tiles = num_tiles_to_flip;

        for (int i = 0; i < num_flipped_tiles; i++)
        {
            last_flipped[i] = tile_idx[i];        
            tile_map[tile_idx[i]] ^= 1;
        }
        compute_P_matrix (tile_map, diff_P);
        // unflip it:
        for (int i = 0; i < num_flipped_tiles; i++)
            tile_map[tile_idx[i]] ^= 1;        

        #pragma omp parallel for
        for (int ii = 0; ii < yee_sz; ii++)
            diff_P[ii] -= yee_grid[ii];

        num_yee_flipped = 0;
        {
        int jj = 0;
        for (int ii = 0; ii < yee_sz; ii++)
        {
            if (diff_P[ii] != 0)
            {
                P_indices[num_yee_flipped] = ii;
                if (diff_P[ii]>0)
                {
                    for (int i = 0; i < num_eobj; i++)
                        RU[num_yee_flipped*num_eobj+i] = R(jj+i);
                } else {
                    for (int i = 0; i < num_eobj; i++)                
                        RU[num_yee_flipped*num_eobj+i] = -R(jj+i);
                }
                num_yee_flipped++;
            }
            jj += num_eobj;
        }
        }

        // obtain V=Q*(G-I) part: (store as V^T column major)
        // this is just a copy taking O(M*N) time
        #pragma omp parallel for
        for (int jj = 0; jj < num_yee_flipped; jj++)
        {
            int base = jj*yee_sz;
            int base_G = P_indices[jj];
            for (int ii = 0; ii < yee_sz; ii++)
            {
                V.data()[base+ii] = G.data()[base_G+ii*yee_sz];
            }
            V.data()[base+P_indices[jj]] -= 1.0_d;
        }

        // obtain V*inv(C)*U+I part:
        #pragma omp parallel for
        for (int ii = 0; ii < num_yee_flipped*num_yee_flipped; ii++)
        {
            int col = ii / num_yee_flipped;
            int row = ii - col * num_yee_flipped;

            int base_V = row*yee_sz;
            int base_C = P_indices[col]*yee_sz;

            std::complex<precision> ctmp;
            ctmp=0;
            for (int i = 0; i < yee_sz; i++)
                ctmp += V.data()[base_V+i]*C.data()[base_C+i];
            if (diff_P[P_indices[col]]<0)
                ctmp = -ctmp;

            if (row == col)
                ctmp += 1.0_d;                
            VinvCU[ii] = ctmp;

        }

        // get: Snew = inv(C)*P_n*Einc = S + inv(C)*dP*Einc
        #pragma omp parallel for
        for (int ii = 0; ii < yee_sz; ii++)
        {
            std::complex<precision> ctmp = 0;
            for (int jj = 0; jj < num_yee_flipped; jj++)
                ctmp += C.data()[ii+P_indices[jj]*yee_sz]*Einc.data()[P_indices[jj]]*(precision)diff_P[P_indices[jj]];
            Snew.data()[ii] = S.data()[ii] - ctmp;
        }
        // now multiply V*Snew:
        #pragma omp parallel for
        for (int ii = 0; ii < num_yee_flipped; ii++)
        {
            int base;

            base = ii*yee_sz;
            std::complex<precision> ctmp;
            ctmp=0;
            for (int jj = 0; jj < yee_sz; jj++)
            {
                ctmp += V.data()[base+jj]*Snew.data()[jj];
            }
            VSnew[ii] = ctmp;
        }

        // now compute R*inv[I+V*inv(C)*U:
        int nrhs = 1; // number of right-hand sides (we have a single vector b)
        int lda = num_yee_flipped;  // leading dimension of A
        int ldb = num_yee_flipped;  // leading dimension of b (or x)
        int info = 0; // output status

#ifdef USE_SINGLE_PRECISION
        cgetrf_(&num_yee_flipped, &num_yee_flipped,
            reinterpret_cast<lp_complex*>(&VinvCU[0]), &num_yee_flipped,
            &ipiv[0], &info);
#else
        zgetrf_(&num_yee_flipped, &num_yee_flipped,
            reinterpret_cast<lp_complex*>(&VinvCU[0]), &num_yee_flipped,
            &ipiv[0], &info);
#endif            
        // Solve the system using the computed LU factors.
        char trans = 'N';  // 'N' indicates no transpose.
#ifdef USE_SINGLE_PRECISION
        cgetrs_(&trans, &num_yee_flipped, &nrhs,
            reinterpret_cast<lp_complex*>(&VinvCU[0]), &num_yee_flipped,
            &ipiv[0],
            reinterpret_cast<lp_complex*>(&VSnew[0]), &num_yee_flipped,
            &info);
#else
        zgetrs_(&trans, &num_yee_flipped, &nrhs,
            reinterpret_cast<lp_complex*>(&VinvCU[0]), &num_yee_flipped,
            &ipiv[0],
            reinterpret_cast<lp_complex*>(&VSnew[0]), &num_yee_flipped,
            &info);
#endif

        
        #pragma omp parallel for
        for (int j = 0; j < num_eobj; j++)
        {
            std::complex<precision> ctmp;
            //#pragma omp parallel for reduction(+:ctmp)
            ctmp = 0;
            for (int i = 0; i < yee_sz; i++)
            {
                ctmp += Gobj(j,i)*Snew.data()[i];
            }
            newEobj(j) = ctmp;
        }

        #pragma omp parallel for
        for (int j = 0; j < num_eobj; j++)
        {
            std::complex<precision> ctmp;        
            ctmp = 0;
            for (int i = 0; i < num_yee_flipped; i++)
            {
                ctmp += RU[i*num_eobj+j]*VSnew[i];
            }
            newEobj(j) -= ctmp;
        }
        return newEobj(0);
    }

    // compute P * Einc from paper:
    void compute_PEinc (void)
    {
        #pragma omp parallel for
        for (int ii = 0; ii < yee_sz; ii++)
        {
            if (yee_grid[ii])
                PEinc.data()[ii] = -Einc.data()[ii];
            else
                PEinc.data()[ii] = 0.0;
        }
    }

    // Select which tile on the grid to be used as the source port.
    // The port_y_offset and port_y_size variables are specified in
    // in terms of the number of Yee cells. The default are offset 1
    // and size 2. Note: Currently, only x directed ports are supported,
    // although adding y directed ports would be straightforward.
    void set_source_port (int port_tile_x, int port_tile_y, int port_y_offset, int port_y_size)
    {
        if (port_y_offset + port_y_size > tiles_yee_y)
        {
            std::cout << "set_source_port: invalid port y position and height parameters." << std::endl;
            exit(0);
        }
        port_height = port_y_size;

        int index;

        // The first row of Gobj is going to be used to compute the voltage across the port by
        // summing the E field values corresponding to the port location:
        Gobj.row(0).setZero();
        // This is the incident E field (Einc) produced by the current source densities comprising
        // the source port, Jsrc.
        Einc.setZero();
        // These are used to compute the component of the two xobj values produced by the incident
        // port source excitation needed to compute the Near-Field to Far-Field Transformation for
        // computing the directivity (D).
        LN1inc = 0;
        LN2inc = 0;
        //std::cout << "got here." << std::endl;
        index=0;
        port_indices.resize(port_width*port_height);
        for (int jj = 0; jj < port_height; jj++)
        {
            for (int ii = 0; ii < port_width; ii++)
            {
                int base_x = (port_tile_x*tiles_yee_x+ii) + (port_tile_y*tiles_yee_y+jj+port_y_offset)*yee_x_xsz;
                port_indices[index] = base_x;

                // Gobj is the sum of all the Ex field components of the port taken from the G matrix: 
                Gobj.row(0) += G.row(base_x);
                // Einc is the sum of the fields produced by all the current source J's comprising the port:
                Einc += G.col(base_x);
                // These are quantities extracted from the NF2FF transformation rows in the Gobj matrix
                // and correspond to the part of xobj produced by the incident sources J which must be summed
                // to the xobj contributed by the solved Jopt current densities later and then used to
                // compute the directivity of the structure:
                if (Gobj.rows()>1)
                {
                    LN1inc += Gobj(1,base_x);
                    LN2inc += Gobj(2,base_x);
                }
                index++;
            }
        }
        // Voltage across the port due to the Einc. This must be summed to the voltage due to Jopt to
        // correctly compute the total voltage across the port:
        Vinc = 0;
        for (int jj = 0; jj < port_height; jj++)
        {
            for (int ii = 0; ii < port_width; ii++)
            {
                int base_x = (port_tile_x*tiles_yee_x+ii) + (port_tile_y*tiles_yee_y+jj+port_y_offset)*yee_x_xsz;
                Vinc += Einc(base_x);
            }
        } 
    }

    // This computes all of P from scratch based on a 0-1 supplied tile mapping.
    // It has a bit of overlap for adjacent tiles which share an edge, but the
    // additional computational cost should be very minimal and likely a worthwhile
    // tradeoff for improved code readability.
    void compute_P_matrix (std::vector<int> &input_map, std::vector<int> &output_P)
    {
        int i;

        #pragma omp parallel for
        for (int index = 0; index < yee_sz; index++)
            output_P[index] = 0;

        #pragma omp parallel for
        for (int index = 0; index < map_sz; index++)
        {
            int tx, ty, base_x, base_y;

            ty = index / map_xsz;
            tx = index - ty*map_xsz;

            base_x = tx*tiles_yee_x + ty*tiles_yee_y*yee_x_xsz;
            base_y = yee_x_sz + tx*tiles_yee_x + ty*tiles_yee_y*yee_y_xsz;

            // check if the tile is set:
            if (input_map[index])
            {
                // Set the X yee cells:
                for (int jj = 0; jj < tiles_yee_y+1; jj++)
                {
                    for (int ii = 0; ii < tiles_yee_x; ii++)
                    {
                        output_P[base_x + ii + jj*yee_x_xsz] = 1;
                    }
                }
                // Set the Y yee cells:
                for (int jj = 0; jj < tiles_yee_y; jj++)
                {
                    for (int ii = 0; ii < tiles_yee_x+1; ii++)
                    {
                        output_P[base_y + ii + jj*yee_y_xsz] = 1;
                    }
                }
            }
       }
    }

    // Compute the Yee cell grid representation of an input tilemap
    // and store it internally. This computes the P matrix from the
    // paper (called "yee_grid" in this code).
    void set_tile_map (std::vector<int> &input_map)
    {
        for (int i = 0; i < map_sz; i++)
            tile_map[i] = input_map[i];
        compute_P_matrix (tile_map, yee_grid);
    }

    ~PNGF_Optimizer()
    {
    }
};