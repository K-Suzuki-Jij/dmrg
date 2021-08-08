//
//  Created by Kohei Suzuki on 2020/12/27.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include <iostream>
#include "DMRG.hpp"

void DMRG_Matrix_Vector_Product(const std::string Mat_Type, const CRS &M, const std::vector<double> &V_In, std::vector<double> &V_Out, const std::vector<int> &Inv, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const int p_threads) {
   
   int dim_LLLRRRRL = (int)Basis_LLLRRRRL.LL.size();
   int base_LRRRRL  = Basis_LLLRRRRL.base_LRRRRL;
   int base_RRRL    = Basis_LLLRRRRL.base_RRRL;
   int base_RL      = Basis_LLLRRRRL.base_RL;
   
   if ((int)V_Out.size() != dim_LLLRRRRL) {
      V_Out.resize(dim_LLLRRRRL);
   }
   
   int c_LL = (Mat_Type == "LL");
   int c_LR = (Mat_Type == "LR");
   int c_RR = (Mat_Type == "RR");
   int c_RL = (Mat_Type == "RL");

#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      int LL = Basis_LLLRRRRL.LL[i];
      int LR = Basis_LLLRRRRL.LR[i];
      int RR = Basis_LLLRRRRL.RR[i];
      int RL = Basis_LLLRRRRL.RL[i];
      double val = 0.0;
      
      if (c_LL) {
         for (long j = M.Row[LL]; j < M.Row[LL+1]; j++) {
            int inv = Inv[(long)M.Col[j]*base_LRRRRL + LR*base_RRRL + RR*base_RL + RL];
            if (inv >= 0) {
               val += V_In[inv]*M.Val[j];
            }
         }
      }
      else if (c_LR) {
         for (long j = M.Row[LR]; j < M.Row[LR+1]; j++) {
            int inv = Inv[(long)LL*base_LRRRRL + M.Col[j]*base_RRRL + RR*base_RL + RL];
            if (inv >= 0) {
               val += V_In[inv]*M.Val[j];
            }
         }
      }
      else if (c_RR) {
         for (long j = M.Row[RR]; j < M.Row[RR+1]; j++) {
            int inv = Inv[(long)LL*base_LRRRRL + LR*base_RRRL + M.Col[j]*base_RL + RL];
            if (inv >= 0) {
               val += V_In[inv]*M.Val[j];
            }
         }
      }
      else if (c_RL) {
         for (long j = M.Row[RL]; j < M.Row[RL+1]; j++) {
            int inv = Inv[(long)LL*base_LRRRRL + LR*base_RRRL + RR*base_RL + M.Col[j]];
            if (inv >= 0) {
               val += V_In[inv]*M.Val[j];
            }
         }
      }
      else {
         std::cout << "Error in DMRG_Matrix_Vector_Product" << std::endl;
         std::exit(0);
      }
      V_Out[i] = val;
   }

}

