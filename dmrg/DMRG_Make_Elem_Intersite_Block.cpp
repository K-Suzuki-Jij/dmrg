//
//  Created by Kohei Suzuki on 2020/12/29.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include <iostream>
#include "DMRG.hpp"

void DMRG_Make_Elem_Intersite_Block(std::string Mat_Type, int basis_change_1, int basis_change_2, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv, const CRS &M_1, const CRS &M_2, double coeff, const std::vector<int> &Ele_LL_or_LR_or_RR, std::string Sign_Flag) {
   
   //pow(10,-15)
   if (std::abs(coeff) < 0.000000000000001) {
      return;
   }
   
   int dim_onsite = A_Basis.base_RL;
   int dim_RR     = A_Basis.dim_RR;
   int elem_num   = A_Basis.elem_num;
   
   int c_LLLR = (Mat_Type  == "LL_LR");
   int c_LLRR = (Mat_Type  == "LL_RR");
   int c_LLRL = (Mat_Type  == "LL_RL");
   int c_LRRR = (Mat_Type  == "LR_RR");
   int c_LRRL = (Mat_Type  == "LR_RL");
   int c_RRRL = (Mat_Type  == "RR_RL");
   
   int sign = 1;
   if (Sign_Flag == "Yes" && Ele_LL_or_LR_or_RR[basis_change_1]%2 == 0) {
      sign = -1;
   }
   
   for (long i = M_1.Row[basis_change_1]; i < M_1.Row[basis_change_1 + 1]; i++) {
      int    col_1 = M_1.Col[i];
      double val_1 = M_1.Val[i];
      for (long j = M_2.Row[basis_change_2]; j < M_2.Row[basis_change_2 + 1]; j++) {
         int inv = -1;
         
         if (c_LLRR || c_LRRR) {
            inv = Inv[col_1*dim_RR + M_2.Col[j]];
         }
         else if (c_LLLR || c_LLRL || c_LRRL || c_RRRL) {
            inv = Inv[col_1*dim_onsite + M_2.Col[j]];
         }
         else {
            std::cout << "Error in DMRG_Make_Elem_Intersite_Block at 1" << std::endl;
            std::exit(0);
         }
         
         if (inv < 0) {
            std::cout << "Error in DMRG_Make_Elem_Intersite_Block at 2" << std::endl;
            std::exit(0);
         }
         
         int check = A_Basis.Check_Basis[inv];
         
         if (check == -1) {
            A_Basis.Check_Basis[inv] = elem_num;
            A_Basis.Val[elem_num]    = sign*coeff*val_1*M_2.Val[j];
            A_Basis.Inv[elem_num]    = inv;
            elem_num++;
         }
         else {
            A_Basis.Val[check] += sign*coeff*val_1*M_2.Val[j];
         }
      }
   }
   
   A_Basis.elem_num = elem_num;
   
}
