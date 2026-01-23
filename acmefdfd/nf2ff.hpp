// Classes for near field to far field transformations
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


// make integration matrices using Simpson's rule
class simpson_integrator
{
    public:
        simpson_integrator(void)
        {
            return;
        }

        // generate a 1D matrix to approximate the integral of a function
        // over the interval [x1,x2] using Simpson's rule
        Eigen::VectorXd get_integration_matrix_1D (precision x1, precision x2, int n_points)
        {
            Eigen::VectorXd integMat;

            if (!(n_points&1))
            {
                std::cout << "Number of intervals must be odd for Simpson's rule!" << std::endl;
                return (Eigen::VectorXd)0;
            }

            integMat.resize(n_points);

            integMat(0) = 1.0_d;
            integMat(n_points-1) = 1.0_d;
            for (int i = 1; i < n_points-1; i++)
            {
                if (i&1)
                    integMat(i) = 4.0_d;
                else
                    integMat(i) = 2.0_d;
            }
            integMat = integMat * (x2-x1)/(3.0_d*(precision)(n_points-1));
            return integMat;
        }

        // generate a 2D matrix to approximate the integral of a function
        // over the rectangle [x1,x2] x [y1,y2] using Simpson's rule
        Eigen::VectorXd get_integration_matrix_2D (precision x1, precision x2, precision y1, precision y2, int n_points_x, int n_points_y)
        {
            Eigen::VectorXd integMat;
            precision weight;
            int index;

            if ((!(n_points_x&1))||(!(n_points_y&1)))
            {
                std::cout << "Number of intervals must be odd for Simpson's rule!" << std::endl;
                return (Eigen::VectorXd)0;
            }

            integMat.resize(n_points_x*n_points_y);

            index=0;
            for (int j = 0; j < n_points_y; j++)
            {
                if ((j==0)||(j==(n_points_y-1)))
                    weight = 1.0_d;
                else if (j&1)
                    weight = 4.0_d;
                else
                    weight = 2.0_d;
                for (int i = 0; i < n_points_x; i++)
                {
                    if ((i==0)||(i==(n_points_x-1)))
                        integMat(index) = 1.0_d*weight;
                    else if (i&1)
                        integMat(index) = 4.0_d*weight;
                    else
                        integMat(index) = 2.0_d*weight;
                    index++;
                }
            }

            integMat = integMat * (x2-x1)*(y2-y1)/(9.0_d*(precision)(n_points_x-1)*(precision)(n_points_y-1));
            return integMat;
        }

        ~simpson_integrator()
        {
            return;
        }        
};


// near field to far field transformations
// Based on Computational Electrodynamics by Taflove and Hagness
// Chapter 08: Near-to-Far-Field Transformation
class NFtoFF
{
    public:
        precision k0;
        int num_theta, num_phi;
        int nfbox_xsz, nfbox_ysz, nfbox_zsz;
        precision step_dx, step_dy, step_dz;
        std::vector<precision> theta_arr, phi_arr;
        std::vector<std::vector<std::complex<precision>>> rcospsi_xP, rcospsi_xN;
        std::vector<std::vector<std::complex<precision>>> rcospsi_yP, rcospsi_yN;
        std::vector<std::vector<std::complex<precision>>> rcospsi_zP, rcospsi_zN;
        simpson_integrator srint;
        Eigen::VectorXd integMat_x, integMat_y, integMat_z;

        NFtoFF(precision theta_start, precision theta_end,
               precision phi_start, precision phi_end,
               int tnum_theta, int tnum_phi,
               precision tstep_dx, precision tstep_dy, precision tstep_dz,
               int tnfbox_xsz, int tnfbox_ysz, int tnfbox_zsz,
               precision tk0)
        {
            int index, index2;
            precision xpos, ypos, zpos;
            precision cur_th, cur_phi;

            // free-space wavenumber:
            k0 = tk0;

            num_theta = tnum_theta;
            num_phi = tnum_phi;

            step_dx = tstep_dx;
            step_dy = tstep_dy;
            step_dz = tstep_dz;

            nfbox_xsz = tnfbox_xsz;
            nfbox_ysz = tnfbox_ysz;
            nfbox_zsz = tnfbox_zsz;

            integMat_x = srint.get_integration_matrix_2D(-(precision)((nfbox_ysz-1)/2) * step_dy,
                                                          (precision)((nfbox_ysz-1)/2) * step_dy,
                                                         -(precision)((nfbox_zsz-1)/2) * step_dz,
                                                          (precision)((nfbox_zsz-1)/2) * step_dz,
                                                          nfbox_ysz, nfbox_zsz);
            integMat_y = srint.get_integration_matrix_2D(-(precision)((nfbox_xsz-1)/2) * step_dx,
                                                          (precision)((nfbox_xsz-1)/2) * step_dx,
                                                         -(precision)((nfbox_zsz-1)/2) * step_dz,
                                                          (precision)((nfbox_zsz-1)/2) * step_dz,
                                                          nfbox_xsz, nfbox_zsz);
            integMat_z = srint.get_integration_matrix_2D(-(precision)((nfbox_xsz-1)/2) * step_dx,
                                                          (precision)((nfbox_xsz-1)/2) * step_dx,
                                                         -(precision)((nfbox_ysz-1)/2) * step_dy,
                                                          (precision)((nfbox_ysz-1)/2) * step_dy,
                                                          nfbox_xsz, nfbox_ysz);
                                                       

            theta_arr.resize(num_theta);
            phi_arr.resize(num_phi);
            rcospsi_xP.resize(num_theta*num_phi);
            rcospsi_xN.resize(num_theta*num_phi);
            rcospsi_yP.resize(num_theta*num_phi);
            rcospsi_yN.resize(num_theta*num_phi);           
            rcospsi_zP.resize(num_theta*num_phi);
            rcospsi_zN.resize(num_theta*num_phi);

            if (num_theta>1)
            {
                for (int i = 0; i < num_theta; i++)
                    theta_arr[i] = theta_start + (theta_end-theta_start)*(precision)i/(precision)(num_theta-1);
            } else
                theta_arr[0] = theta_start;

            if (num_phi>1)
            {
                for (int i = 0; i < num_phi; i++)
                    phi_arr[i] = phi_start + (phi_end-phi_start)*(precision)i/(precision)(num_phi-1);
            } else
                phi_arr[0] = phi_start;                

            // precompute r' cos psi:
            index=0;
            for (int phi_idx = 0; phi_idx < num_phi; phi_idx++)
            {
                cur_phi = phi_arr[phi_idx];
                for (int theta_idx = 0; theta_idx < num_theta; theta_idx++)
                {
                    cur_th = theta_arr[theta_idx];
                    //std::cout << "(th,ph): " << cur_th*(180/M_PI) << ", " << cur_phi*(180/M_PI) << std::endl;

                    rcospsi_xP[index].resize(nfbox_ysz*nfbox_zsz);
                    rcospsi_xN[index].resize(nfbox_ysz*nfbox_zsz);
                    rcospsi_yP[index].resize(nfbox_xsz*nfbox_zsz);
                    rcospsi_yN[index].resize(nfbox_xsz*nfbox_zsz);
                    rcospsi_zP[index].resize(nfbox_xsz*nfbox_ysz);
                    rcospsi_zN[index].resize(nfbox_xsz*nfbox_ysz);

                    xpos = (precision)((nfbox_xsz-1)/2) * step_dx;
                    index2=0;
                    for (int jj = 0; jj < nfbox_zsz; jj++)
                    {
                        zpos = (precision)(-(nfbox_zsz-1)/2+jj) * step_dz;
                        for (int ii = 0; ii < nfbox_ysz; ii++)
                        {
                            ypos = (precision)(-(nfbox_ysz-1)/2+ii) * step_dy;
                            rcospsi_xN[index][index2] = exp(I1*k0*((-xpos)*sin(cur_th)*cos(cur_phi)+
                                                                     ypos*sin(cur_th)*sin(cur_phi)+
                                                                     zpos*cos(cur_th)));
                            rcospsi_xP[index][index2] = exp(I1*k0*(xpos*sin(cur_th)*cos(cur_phi)+
                                                                   ypos*sin(cur_th)*sin(cur_phi)+
                                                                   zpos*cos(cur_th)));                                                                   
                            index2++;
                        }
                    }
                    ypos = (precision)((nfbox_ysz-1)/2) * step_dy;
                    index2=0;
                    for (int jj = 0; jj < nfbox_zsz; jj++)
                    {
                        zpos = (precision)(-(nfbox_zsz-1)/2+jj) * step_dz;
                        for (int ii = 0; ii < nfbox_xsz; ii++)
                        {
                            xpos = (precision)(-(nfbox_xsz-1)/2+ii) * step_dx;
                            rcospsi_yN[index][index2] = exp(I1*k0*(xpos*sin(cur_th)*cos(cur_phi)+
                                                                 (-ypos)*sin(cur_th)*sin(cur_phi)+
                                                                   zpos*cos(cur_th)));
                            rcospsi_yP[index][index2] = exp(I1*k0*(xpos*sin(cur_th)*cos(cur_phi)+
                                                                   ypos*sin(cur_th)*sin(cur_phi)+
                                                                   zpos*cos(cur_th)));                                                                   
                            index2++;
                        }
                    }   
                    zpos = (precision)((nfbox_zsz-1)/2) * step_dz;
                    index2=0;
                    for (int jj = 0; jj < nfbox_ysz; jj++)
                    {
                        ypos = (precision)(-(nfbox_ysz-1)/2+jj) * step_dy;
                        for (int ii = 0; ii < nfbox_xsz; ii++)
                        {
                            xpos = (precision)(-(nfbox_xsz-1)/2+ii) * step_dx;
                            rcospsi_zN[index][index2] = exp(I1*k0*(xpos*sin(cur_th)*cos(cur_phi)+
                                                                   ypos*sin(cur_th)*sin(cur_phi)+
                                                                 (-zpos)*cos(cur_th)));
                            rcospsi_zP[index][index2] = exp(I1*k0*(xpos*sin(cur_th)*cos(cur_phi)+
                                                                   ypos*sin(cur_th)*sin(cur_phi)+
                                                                   zpos*cos(cur_th)));                                                                   
                            index2++;
                        }
                    }                    
                    index++;
                }
            }
        }

        // compute Ntheta, Nphi, Ltheta, Lphi as described in Chapter 8 of Taflove/Hagness:
        void computeRCS (Eigen::VectorXcd &Ey_xP, Eigen::VectorXcd &Ez_xP, Eigen::VectorXcd &Hy_xP, Eigen::VectorXcd &Hz_xP,
                        Eigen::VectorXcd &Ey_xN, Eigen::VectorXcd &Ez_xN, Eigen::VectorXcd &Hy_xN, Eigen::VectorXcd &Hz_xN,
                        Eigen::VectorXcd &Ex_yP, Eigen::VectorXcd &Ez_yP, Eigen::VectorXcd &Hx_yP, Eigen::VectorXcd &Hz_yP,
                        Eigen::VectorXcd &Ex_yN, Eigen::VectorXcd &Ez_yN, Eigen::VectorXcd &Hx_yN, Eigen::VectorXcd &Hz_yN,
                        Eigen::VectorXcd &Ex_zP, Eigen::VectorXcd &Ey_zP, Eigen::VectorXcd &Hx_zP, Eigen::VectorXcd &Hy_zP,
                        Eigen::VectorXcd &Ex_zN, Eigen::VectorXcd &Ey_zN, Eigen::VectorXcd &Hx_zN, Eigen::VectorXcd &Hy_zN,                        
                        Eigen::VectorXd &rcs)
        {
            int index, index2;
            precision cur_th, cur_phi;

            Eigen::VectorXcd Nth_x(Ey_xP.rows()), Nph_x(Ey_xP.rows()), Lth_x(Ey_xP.rows()), Lph_x(Ey_xP.rows());
            Eigen::VectorXcd Nth_y(Ex_yP.rows()), Nph_y(Ex_yP.rows()), Lth_y(Ex_yP.rows()), Lph_y(Ex_yP.rows());
            Eigen::VectorXcd Nth_z(Ex_zP.rows()), Nph_z(Ex_zP.rows()), Lth_z(Ex_zP.rows()), Lph_z(Ex_zP.rows());

            std::complex<precision> int_Nth_xN, int_Nph_xN, int_Lth_xN, int_Lph_xN;
            std::complex<precision> int_Nth_xP, int_Nph_xP, int_Lth_xP, int_Lph_xP;
            std::complex<precision> int_Nth_yN, int_Nph_yN, int_Lth_yN, int_Lph_yN;
            std::complex<precision> int_Nth_yP, int_Nph_yP, int_Lth_yP, int_Lph_yP;
            std::complex<precision> int_Nth_zN, int_Nph_zN, int_Lth_zN, int_Lph_zN;
            std::complex<precision> int_Nth_zP, int_Nph_zP, int_Lth_zP, int_Lph_zP;
            std::complex<precision> int_Nth, int_Nph, int_Lth, int_Lph;
            std::complex<precision> Jx, Jy, Jz, Mx, My, Mz;
            std::cout << "num_theta: " << num_theta << ", num_phi: " << num_phi << std::endl;
            index=0;
            for (int phi_idx = 0; phi_idx < num_phi; phi_idx++)
            {
                cur_phi = phi_arr[phi_idx];
                for (int theta_idx = 0; theta_idx < num_theta; theta_idx++)
                {
                    cur_th = theta_arr[theta_idx];
                    
                    //xpos = (precision)((nfbox_xsz-1)/2) * step_dx;
                    index2=0;
                    for (int jj = 0; jj < nfbox_zsz; jj++)
                    {
                        //zpos = (precision)(-(nfbox_zsz-1)/2+jj) * step_dz;
                        for (int ii = 0; ii < nfbox_ysz; ii++)
                        {
                            //ypos = (precision)(-(nfbox_ysz-1)/2+ii) * step_dy;
                            Jx = 0;
                            Jy = Hz_xN(index2);
                            Jz = -Hy_xN(index2);
                            Nth_x(index2) = (Jx*cos(cur_th)*cos(cur_phi)+
                                            Jy*cos(cur_th)*sin(cur_phi)-
                                            Jz*sin(cur_th))*rcospsi_xN[index][index2];
                            Nph_x(index2) = (-Jx*sin(cur_phi)+Jy*cos(cur_phi))*rcospsi_xN[index][index2];
                            Mx = 0;
                            My = -Ez_xN(index2);
                            Mz = Ey_xN(index2);

                            Lth_x(index2) = (Mx*cos(cur_th)*cos(cur_phi)+
                                             My*cos(cur_th)*sin(cur_phi)-
                                             Mz*sin(cur_th))*rcospsi_xN[index][index2];
                            Lph_x(index2) = (-Mx*sin(cur_phi)+My*cos(cur_phi))*rcospsi_xN[index][index2];
                            index2++;
                        }
                    }
                    int_Nth_xN = integMat_x.transpose()*Nth_x;
                    int_Nph_xN = integMat_x.transpose()*Nph_x;
                    int_Lth_xN = integMat_x.transpose()*Lth_x;
                    int_Lph_xN = integMat_x.transpose()*Lph_x;
                    index2=0;
                    for (int jj = 0; jj < nfbox_zsz; jj++)
                    {
                        //zpos = (precision)(-(nfbox_zsz-1)/2+jj) * step_dz;
                        for (int ii = 0; ii < nfbox_ysz; ii++)
                        {
                            //ypos = (precision)(-(nfbox_ysz-1)/2+ii) * step_dy;
                            Jx = 0;
                            Jy = -Hz_xP(index2);
                            Jz = Hy_xP(index2);
                            Nth_x(index2) = (Jx*cos(cur_th)*cos(cur_phi)+
                                            Jy*cos(cur_th)*sin(cur_phi)-
                                            Jz*sin(cur_th))*rcospsi_xP[index][index2];
                            Nph_x(index2) = (-Jx*sin(cur_phi)+Jy*cos(cur_phi))*rcospsi_xP[index][index2];
                            Mx = 0;
                            My = Ez_xP(index2);
                            Mz = -Ey_xP(index2);

                            Lth_x(index2) = (Mx*cos(cur_th)*cos(cur_phi)+
                                             My*cos(cur_th)*sin(cur_phi)-
                                             Mz*sin(cur_th))*rcospsi_xP[index][index2];
                            Lph_x(index2) = (-Mx*sin(cur_phi)+My*cos(cur_phi))*rcospsi_xP[index][index2];
                            index2++;
                        }
                    }
                    int_Nth_xP = integMat_x.transpose()*Nth_x;
                    int_Nph_xP = integMat_x.transpose()*Nph_x;
                    int_Lth_xP = integMat_x.transpose()*Lth_x;
                    int_Lph_xP = integMat_x.transpose()*Lph_x;
                    index2=0;
                    for (int jj = 0; jj < nfbox_zsz; jj++)
                    {
                        //zpos = (precision)(-(nfbox_zsz-1)/2+jj) * step_dz;
                        for (int ii = 0; ii < nfbox_xsz; ii++)
                        {
                            //xpos = (precision)(-(nfbox_xsz-1)/2+ii) * step_dx;
                            Jx = -Hz_yN(index2);
                            Jy = 0;
                            Jz = Hx_yN(index2);
                            Nth_y(index2) = (Jx*cos(cur_th)*cos(cur_phi)+
                                            Jy*cos(cur_th)*sin(cur_phi)-
                                            Jz*sin(cur_th))*rcospsi_yN[index][index2];
                            Nph_y(index2) = (-Jx*sin(cur_phi)+Jy*cos(cur_phi))*rcospsi_yN[index][index2];
                            Mx = Ez_yN(index2);
                            My = 0;
                            Mz = -Ex_yN(index2);

                            Lth_y(index2) = (Mx*cos(cur_th)*cos(cur_phi)+
                                             My*cos(cur_th)*sin(cur_phi)-
                                             Mz*sin(cur_th))*rcospsi_yN[index][index2];
                            Lph_y(index2) = (-Mx*sin(cur_phi)+My*cos(cur_phi))*rcospsi_yN[index][index2];
                            index2++;
                        }
                    }
                    int_Nth_yN = integMat_y.transpose()*Nth_y;
                    int_Nph_yN = integMat_y.transpose()*Nph_y;
                    int_Lth_yN = integMat_y.transpose()*Lth_y;
                    int_Lph_yN = integMat_y.transpose()*Lph_y;     
                    index2=0;
                    for (int jj = 0; jj < nfbox_zsz; jj++)
                    {
                        //zpos = (precision)(-(nfbox_zsz-1)/2+jj) * step_dz;
                        for (int ii = 0; ii < nfbox_xsz; ii++)
                        {
                            //xpos = (precision)(-(nfbox_xsz-1)/2+ii) * step_dx;
                            Jx = Hz_yP(index2);
                            Jy = 0;
                            Jz = -Hx_yP(index2);
                            Nth_y(index2) = (Jx*cos(cur_th)*cos(cur_phi)+
                                            Jy*cos(cur_th)*sin(cur_phi)-
                                            Jz*sin(cur_th))*rcospsi_yP[index][index2];
                            Nph_y(index2) = (-Jx*sin(cur_phi)+Jy*cos(cur_phi))*rcospsi_yP[index][index2];
                            Mx = -Ez_yP(index2);
                            My = 0;
                            Mz = Ex_yP(index2);

                            Lth_y(index2) = (Mx*cos(cur_th)*cos(cur_phi)+
                                             My*cos(cur_th)*sin(cur_phi)-
                                             Mz*sin(cur_th))*rcospsi_yP[index][index2];
                            Lph_y(index2) = (-Mx*sin(cur_phi)+My*cos(cur_phi))*rcospsi_yP[index][index2];
                            index2++;
                        }
                    }
                    int_Nth_yP = integMat_y.transpose()*Nth_y;
                    int_Nph_yP = integMat_y.transpose()*Nph_y;
                    int_Lth_yP = integMat_y.transpose()*Lth_y;
                    int_Lph_yP = integMat_y.transpose()*Lph_y;   
                    index2=0;
                    for (int jj = 0; jj < nfbox_ysz; jj++)
                    {
                        //ypos = (precision)(-(nfbox_ysz-1)/2+jj) * step_dy;
                        for (int ii = 0; ii < nfbox_xsz; ii++)
                        {
                            //xpos = (precision)(-(nfbox_xsz-1)/2+ii) * step_dx;
                            Jx = Hy_zN(index2);
                            Jy = -Hx_zN(index2);
                            Jz = 0;
                            Nth_z(index2) = (Jx*cos(cur_th)*cos(cur_phi)+
                                            Jy*cos(cur_th)*sin(cur_phi)-
                                            Jz*sin(cur_th))*rcospsi_zN[index][index2];
                            Nph_z(index2) = (-Jx*sin(cur_phi)+Jy*cos(cur_phi))*rcospsi_zN[index][index2];
                            Mx = -Ey_zN(index2);
                            My = Ex_zN(index2);
                            Mz = 0;
                            Lth_z(index2) = (Mx*cos(cur_th)*cos(cur_phi)+
                                             My*cos(cur_th)*sin(cur_phi)-
                                             Mz*sin(cur_th))*rcospsi_zN[index][index2];
                            Lph_z(index2) = (-Mx*sin(cur_phi)+My*cos(cur_phi))*rcospsi_zN[index][index2];
                            index2++;
                        }
                    }
                    int_Nth_zN = integMat_z.transpose()*Nth_z;
                    int_Nph_zN = integMat_z.transpose()*Nph_z;
                    int_Lth_zN = integMat_z.transpose()*Lth_z;
                    int_Lph_zN = integMat_z.transpose()*Lph_z;                      
                    index2=0;
                    for (int jj = 0; jj < nfbox_ysz; jj++)
                    {
                        //ypos = (precision)(-(nfbox_ysz-1)/2+jj) * step_dy;
                        for (int ii = 0; ii < nfbox_xsz; ii++)
                        {
                            //xpos = (precision)(-(nfbox_xsz-1)/2+ii) * step_dx;
                            Jx = -Hy_zP(index2);
                            Jy = Hx_zP(index2);
                            Jz = 0;
                            Nth_z(index2) = (Jx*cos(cur_th)*cos(cur_phi)+
                                            Jy*cos(cur_th)*sin(cur_phi)-
                                            Jz*sin(cur_th))*rcospsi_zP[index][index2];
                            Nph_z(index2) = (-Jx*sin(cur_phi)+Jy*cos(cur_phi))*rcospsi_zP[index][index2];
                            Mx = Ey_zP(index2);
                            My = -Ex_zP(index2);
                            Mz = 0;
                            Lth_z(index2) = (Mx*cos(cur_th)*cos(cur_phi)+
                                             My*cos(cur_th)*sin(cur_phi)-
                                             Mz*sin(cur_th))*rcospsi_zP[index][index2];
                            Lph_z(index2) = (-Mx*sin(cur_phi)+My*cos(cur_phi))*rcospsi_zP[index][index2];
                            index2++;
                        }
                    }
                    int_Nth_zP = integMat_z.transpose()*Nth_z;
                    int_Nph_zP = integMat_z.transpose()*Nph_z;
                    int_Lth_zP = integMat_z.transpose()*Lth_z;
                    int_Lph_zP = integMat_z.transpose()*Lph_z;
                    
                    int_Nth = int_Nth_xP + int_Nth_xN + int_Nth_yP + int_Nth_yN + int_Nth_zP + int_Nth_zN;
                    int_Nph = int_Nph_xP + int_Nph_xN + int_Nph_yP + int_Nph_yN + int_Nph_zP + int_Nph_zN;
                    int_Lth = int_Lth_xP + int_Lth_xN + int_Lth_yP + int_Lth_yN + int_Lth_zP + int_Lth_zN;
                    int_Lph = int_Lph_xP + int_Lph_xN + int_Lph_yP + int_Lph_yN + int_Lph_zP + int_Lph_zN;

                    precision tmp1, tmp2;
                    tmp1 = std::real((int_Lph+eta0*int_Nth)*std::conj(int_Lph+eta0*int_Nth));
                    tmp2 = std::real((int_Lth-eta0*int_Nph)*std::conj(int_Lth-eta0*int_Nph));

                    rcs(index) = (k0*k0/(32*M_PI*M_PI*eta0))*(tmp1+tmp2);
                    //std::cout << "RCS: " << rcs(index) << std::endl;
                    index++;
                }
            }

        }

        // compute Ntheta, Nphi, Ltheta, Lphi as described in Chapter 8 of Taflove/Hagness:
        void generate_RCS_matrix (Eigen::SparseMatrix<std::complex<precision>> &rcs_mat)
        {
            int index, index2;
            precision cur_th, cur_phi;

            rcs_mat.resize(num_theta*num_phi*2, 4*2*(nfbox_ysz*nfbox_zsz+nfbox_xsz*nfbox_zsz+nfbox_xsz*nfbox_ysz));
            std::vector<Eigen::Triplet<std::complex<precision>>> triplets;
            triplets.reserve(8*(nfbox_ysz*nfbox_zsz*2+nfbox_xsz*nfbox_zsz*2+nfbox_xsz*nfbox_ysz*2)); 
            
            int base_Ey_xP, base_Ez_xP, base_Hy_xP, base_Hz_xP;
            int base_Ey_xN, base_Ez_xN, base_Hy_xN, base_Hz_xN;
            int base_Ex_yP, base_Ez_yP, base_Hx_yP, base_Hz_yP;
            int base_Ex_yN, base_Ez_yN, base_Hx_yN, base_Hz_yN;
            int base_Ex_zP, base_Ey_zP, base_Hx_zP, base_Hy_zP;
            int base_Ex_zN, base_Ey_zN, base_Hx_zN, base_Hy_zN;

            base_Ey_xP = 0;
            base_Ez_xP = base_Ey_xP + nfbox_ysz*nfbox_zsz;
            base_Hy_xP = base_Ez_xP + nfbox_ysz*nfbox_zsz;
            base_Hz_xP = base_Hy_xP + nfbox_ysz*nfbox_zsz;

            base_Ey_xN = base_Hz_xP + nfbox_ysz*nfbox_zsz;
            base_Ez_xN = base_Ey_xN + nfbox_ysz*nfbox_zsz;
            base_Hy_xN = base_Ez_xN + nfbox_ysz*nfbox_zsz;
            base_Hz_xN = base_Hy_xN + nfbox_ysz*nfbox_zsz;

            base_Ex_yP = base_Hz_xN + nfbox_ysz*nfbox_zsz;
            base_Ez_yP = base_Ex_yP + nfbox_xsz*nfbox_zsz;
            base_Hx_yP = base_Ez_yP + nfbox_xsz*nfbox_zsz;
            base_Hz_yP = base_Hx_yP + nfbox_xsz*nfbox_zsz;

            base_Ex_yN = base_Hz_yP + nfbox_xsz*nfbox_zsz;
            base_Ez_yN = base_Ex_yN + nfbox_xsz*nfbox_zsz;
            base_Hx_yN = base_Ez_yN + nfbox_xsz*nfbox_zsz;
            base_Hz_yN = base_Hx_yN + nfbox_xsz*nfbox_zsz;

            base_Ex_zP = base_Hz_yN + nfbox_xsz*nfbox_zsz;
            base_Ey_zP = base_Ex_zP + nfbox_xsz*nfbox_ysz;
            base_Hx_zP = base_Ey_zP + nfbox_xsz*nfbox_ysz;
            base_Hy_zP = base_Hx_zP + nfbox_xsz*nfbox_ysz;

            base_Ex_zN = base_Hy_zP + nfbox_xsz*nfbox_ysz;
            base_Ey_zN = base_Ex_zN + nfbox_xsz*nfbox_ysz;
            base_Hx_zN = base_Ey_zN + nfbox_xsz*nfbox_ysz;
            base_Hy_zN = base_Hx_zN + nfbox_xsz*nfbox_ysz;

            //std::cout << "base end: " << base_Hy_zN + nfbox_xsz*nfbox_ysz << std::endl;
            //std::cout << "num_theta: " << num_theta << ", num_phi: " << num_phi << std::endl;
            index=0;
            for (int phi_idx = 0; phi_idx < num_phi; phi_idx++)
            {
                cur_phi = phi_arr[phi_idx];
                precision cos_phi, sin_phi;

                cos_phi = cos(cur_phi);
                sin_phi = sin(cur_phi);
                for (int theta_idx = 0; theta_idx < num_theta; theta_idx++)
                {
                    cur_th = theta_arr[theta_idx];
                    precision cos_th, sin_th;
                    cos_th = cos(cur_th);
                    sin_th = sin(cur_th);
                    
                    //xpos = (precision)((nfbox_xsz-1)/2) * step_dx;
                    index2=0;
                    for (int jj = 0; jj < nfbox_zsz; jj++)
                    {
                        //zpos = (precision)(-(nfbox_zsz-1)/2+jj) * step_dz;
                        for (int ii = 0; ii < nfbox_ysz; ii++)
                        {
                            //ypos = (precision)(-(nfbox_ysz-1)/2+ii) * step_dy;

                            // for -X:
                            // Nth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hz_xN+index2, 
                                eta0*cos_th*sin_phi*rcospsi_xN[index][index2]*integMat_x(index2)));
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hy_xN+index2, 
                                -eta0*(-sin_th)*rcospsi_xN[index][index2]*integMat_x(index2)));
                            // Lphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Ez_xN+index2, 
                                -cos_phi*rcospsi_xN[index][index2]*integMat_x(index2)));                            

                            // Nphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Hz_xN+index2, 
                                -eta0*cos_phi*rcospsi_xN[index][index2]*integMat_x(index2)));                            
                            // Lth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ez_xN+index2, 
                                -cos_th*sin_phi*rcospsi_xN[index][index2]*integMat_x(index2)));                      
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ey_xN+index2, 
                                (-sin_th)*rcospsi_xN[index][index2]*integMat_x(index2)));

                            // for +X:
                            // Nth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hz_xP+index2, 
                                -eta0*cos_th*sin_phi*rcospsi_xP[index][index2]*integMat_x(index2)));
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hy_xP+index2, 
                                eta0*(-sin_th)*rcospsi_xP[index][index2]*integMat_x(index2)));
                            // Lphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Ez_xP+index2, 
                                cos_phi*rcospsi_xP[index][index2]*integMat_x(index2)));                            

                            // Nphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Hz_xP+index2, 
                                -(-eta0)*cos_phi*rcospsi_xP[index][index2]*integMat_x(index2)));                            
                            // Lth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ez_xP+index2, 
                                cos_th*sin_phi*rcospsi_xP[index][index2]*integMat_x(index2)));                      
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ey_xP+index2, 
                                -(-sin_th)*rcospsi_xP[index][index2]*integMat_x(index2)));

                            index2++;                        
                        }
                    }
                    
                    index2=0;
                    for (int jj = 0; jj < nfbox_zsz; jj++)
                    {
                        //zpos = (precision)(-(nfbox_zsz-1)/2+jj) * step_dz;
                        for (int ii = 0; ii < nfbox_xsz; ii++)
                        {

                            // for -Y:
                            // Nth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hz_yN+index2, 
                                -eta0*cos_th*cos_phi*rcospsi_yN[index][index2]*integMat_y(index2)));
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hx_yN+index2, 
                                eta0*(-sin_th)*rcospsi_yN[index][index2]*integMat_y(index2)));
                            // Lphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Ez_yN+index2, 
                                (-sin_phi)*rcospsi_yN[index][index2]*integMat_y(index2)));                            

                            // Nphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Hz_yN+index2, 
                                -(-eta0)*(-sin_phi)*rcospsi_yN[index][index2]*integMat_y(index2)));                            
                            // Lth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ez_yN+index2, 
                                cos_th*cos_phi*rcospsi_yN[index][index2]*integMat_y(index2)));                      
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ex_yN+index2, 
                                -(-sin_th)*rcospsi_yN[index][index2]*integMat_y(index2)));

                            // for +Y:
                            // Nth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hz_yP+index2, 
                                -(-eta0)*cos_th*cos_phi*rcospsi_yP[index][index2]*integMat_y(index2)));
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hx_yP+index2, 
                                -eta0*(-sin_th)*rcospsi_yP[index][index2]*integMat_y(index2)));
                            // Lphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Ez_yP+index2, 
                                -(-sin_phi)*rcospsi_yP[index][index2]*integMat_y(index2)));                            

                            // Nphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Hz_yP+index2, 
                                -(-(-eta0))*(-sin_phi)*rcospsi_yP[index][index2]*integMat_y(index2)));                            
                            // Lth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ez_yP+index2, 
                                -cos_th*cos_phi*rcospsi_yP[index][index2]*integMat_y(index2)));                      
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ex_yP+index2, 
                                -(-(-sin_th))*rcospsi_yP[index][index2]*integMat_y(index2)));                                

                            index2++;
                        }
                    }

                    index2=0;
                    for (int jj = 0; jj < nfbox_ysz; jj++)
                    {
                        //ypos = (precision)(-(nfbox_ysz-1)/2+jj) * step_dy;
                        for (int ii = 0; ii < nfbox_xsz; ii++)
                        {
                            // for -Z:
                            // Nth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hy_zN+index2, 
                                eta0*cos_th*cos_phi*rcospsi_zN[index][index2]*integMat_z(index2)));
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hx_zN+index2, 
                                -eta0*(cos_th*sin_th)*rcospsi_zN[index][index2]*integMat_z(index2)));
                            // Lphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Ey_zN+index2, 
                                -(-sin_phi)*rcospsi_zN[index][index2]*integMat_z(index2)));                            
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Ex_zN+index2, 
                                (cos_phi)*rcospsi_zN[index][index2]*integMat_z(index2)));

                            // Nphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Hy_zN+index2, 
                                (-eta0)*(-sin_phi)*rcospsi_zN[index][index2]*integMat_z(index2)));  
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Hx_zN+index2, 
                                -(-eta0)*(cos_phi)*rcospsi_zN[index][index2]*integMat_z(index2)));                                                            
                            // Lth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ey_zN+index2, 
                                -cos_th*cos_phi*rcospsi_zN[index][index2]*integMat_z(index2)));                      
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ex_zN+index2, 
                                cos_th*sin_phi*rcospsi_zN[index][index2]*integMat_z(index2)));

                            // for +Z:
                            // Nth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hy_zP+index2, 
                                -eta0*cos_th*cos_phi*rcospsi_zP[index][index2]*integMat_z(index2)));
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Hx_zP+index2, 
                                -(-eta0)*(cos_th*sin_th)*rcospsi_zP[index][index2]*integMat_z(index2)));
                            // Lphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Ey_zP+index2, 
                                -(-(-sin_phi))*rcospsi_zP[index][index2]*integMat_z(index2)));                            
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+0, base_Ex_zP+index2, 
                                -(cos_phi)*rcospsi_zP[index][index2]*integMat_z(index2)));

                            // Nphi:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Hy_zP+index2, 
                                -(-eta0)*(-sin_phi)*rcospsi_zP[index][index2]*integMat_z(index2)));  
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Hx_zP+index2, 
                                -(-(-eta0))*(cos_phi)*rcospsi_zP[index][index2]*integMat_z(index2)));                                                            
                            // Lth:
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ey_zP+index2, 
                                -(-cos_th*cos_phi)*rcospsi_zP[index][index2]*integMat_z(index2)));                      
                            triplets.push_back(Eigen::Triplet<std::complex<precision>>(index*2+1, base_Ex_zP+index2, 
                                -cos_th*sin_phi*rcospsi_zP[index][index2]*integMat_z(index2)));     

                            index2++;
                        }
                    }
                    index++;
                }
            }
                     
            rcs_mat.setFromTriplets(triplets.begin(), triplets.end());
        }

        ~NFtoFF()
        {

        }
};
