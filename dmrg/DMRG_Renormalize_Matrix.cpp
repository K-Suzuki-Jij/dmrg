//
//  Created by Kohei Suzuki on 2021/01/08.
//

#include <iostream>
#include "DMRG.hpp"

void DMRG_Renormalize_Matrix(std::string Mat_Type, const CRS &M, const CCS &Trans_Mat, CRS &M_Out_LL, const DMRG_Block_Information &Block, std::string Sign_Flag, const DMRG_Basis_LLLR &Basis_LLLR, const std::vector<int> &Ele_LL, int p_threads) {
   
   double zero    = 0.000000000000001;//pow(10,-15)
   int dim_renorm = Trans_Mat.col_dim;
   int dim_LLLR   = Trans_Mat.row_dim;
   int dim_LL     = Block.dim_LL;
   int dim_onsite = Block.dim_onsite;
   int c_LL       = Mat_Type == "LL";
   int c_LR       = Mat_Type == "LR";
   int c_sign     = Sign_Flag == "Yes";
   
   if (dim_LL*dim_onsite < dim_LLLR) {
      std::cout << "Error in DMRG_Renormalize_Matrix at 1" << std::endl;
      std::cout << "dim_LL*dim_onsite=" << dim_LL*dim_onsite << ", dim_LLLR=" << dim_LLLR << std::endl;
      std::exit(0);
   }
   
   std::vector<double> Temp_V1(dim_LL*dim_onsite, 0.0);
   std::vector<double> Temp_V2(dim_LLLR         , 0.0);
   CCS Temp_M;
      
   Temp_M.Col.push_back(0);
   
   for (int LL_new = 0; LL_new < dim_renorm; LL_new++) {
      
#pragma omp parallel for num_threads (p_threads)
      for (long i = Trans_Mat.Col[LL_new]; i < Trans_Mat.Col[LL_new+1]; i++) {
         int LL = Basis_LLLR.LL[Trans_Mat.Row[i]];
         int LR = Basis_LLLR.LR[Trans_Mat.Row[i]];
         Temp_V1[LL*dim_onsite + LR] = Trans_Mat.Val[i];
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int LLLR = 0; LLLR < dim_LLLR; LLLR++) {
         int LL = Basis_LLLR.LL[LLLR];
         int LR = Basis_LLLR.LR[LLLR];
         double val = 0.0;
         if (c_LL) {
            for (long i = M.Row[LL]; i < M.Row[LL+1]; i++) {
               val += Temp_V1[M.Col[i]*dim_onsite + LR]*M.Val[i];
            }
         }
         else if (c_LR) {
            int sign = 1;
            if (c_sign && Ele_LL[LL]%2 == 0) {
               sign = -1;
            }
            for (long i = M.Row[LR]; i < M.Row[LR+1]; i++) {
               val += Temp_V1[LL*dim_onsite + M.Col[i]]*M.Val[i]*sign;
            }
         }
         else {
            std::cout << "Error in DMRG_Renormalize_Matrix at 2" << std::endl;
            std::exit(1);
         }
         
         Temp_V2[LLLR] = val;
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int LLLR = 0; LLLR < dim_LLLR; LLLR++) {
         if (std::abs(Temp_V2[LLLR]) > zero) {
#pragma omp critical
            {
               Temp_M.Val.push_back(Temp_V2[LLLR]);
               Temp_M.Row.push_back(LLLR);
            }
         }
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (long i = Trans_Mat.Col[LL_new]; i < Trans_Mat.Col[LL_new+1]; i++) {
         int LL = Basis_LLLR.LL[Trans_Mat.Row[i]];
         int LR = Basis_LLLR.LR[Trans_Mat.Row[i]];
         Temp_V1[LL*dim_onsite + LR] = 0.0;
      }
      
      Temp_M.Col.push_back(Temp_M.Row.size());
      
   }
   
   Temp_M.row_dim = dim_LLLR;
   Temp_M.col_dim = dim_renorm;
   
   Clear_Matrix(M_Out_LL);
   
   M_Out_LL.Row.push_back(0);
   for (int LL_new_row = 0; LL_new_row < dim_renorm; LL_new_row++) {
      
#pragma omp parallel for num_threads (p_threads)
      for (long i = Trans_Mat.Col[LL_new_row]; i < Trans_Mat.Col[LL_new_row+1]; i++) {
         Temp_V1[Trans_Mat.Row[i]] = Trans_Mat.Val[i];
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int LL_new_col = 0; LL_new_col < dim_renorm; LL_new_col++) {
         double val = 0.0;
         for (long i = Temp_M.Col[LL_new_col]; i < Temp_M.Col[LL_new_col + 1]; i++) {
            val += Temp_M.Val[i]*Temp_V1[Temp_M.Row[i]];
         }
         Temp_V2[LL_new_col] = val;
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int LL_new_col = 0; LL_new_col < dim_renorm; LL_new_col++) {
         if (std::abs(Temp_V2[LL_new_col]) > zero) {
#pragma omp critical
            {
               M_Out_LL.Val.push_back(Temp_V2[LL_new_col]);
               M_Out_LL.Col.push_back(LL_new_col);
            }
         }
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (long i = Trans_Mat.Col[LL_new_row]; i < Trans_Mat.Col[LL_new_row+1]; i++) {
         Temp_V1[Trans_Mat.Row[i]] = 0.0;
      }
      
      M_Out_LL.Row.push_back(M_Out_LL.Col.size());
   }
   
   M_Out_LL.row_dim = dim_renorm;
   M_Out_LL.col_dim = dim_renorm;
   
   Sort_Col_CRS(M_Out_LL, p_threads);
   
}
