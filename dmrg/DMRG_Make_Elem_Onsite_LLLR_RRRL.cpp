//
//  Created by Kohei Suzuki on 2020/12/28.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include <iostream>
#include "DMRG.hpp"

void DMRG_Make_Elem_Onsite_LLLR_RRRL(const std::string Mat_Type, int basis_change, int basis_no_change, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv, const CRS &M, double coeff) {
   
   //pow(10,-15)
   if (std::abs(coeff) < 0.000000000000001) {
      return;
   }
   
   int elem_num   = A_Basis.elem_num;
   int dim_onsite = A_Basis.base_RL;
   
   int c_1 = (Mat_Type == "LL" || Mat_Type == "RR");
   int c_2 = (Mat_Type == "LR" || Mat_Type == "RL");
   
   for (long i = M.Row[basis_change]; i < M.Row[basis_change + 1]; i++) {
   
      int inv = -1;
      if (c_1) {
         inv = Inv[M.Col[i]*dim_onsite + basis_no_change];
      }
      else if (c_2) {
         inv = Inv[basis_no_change*dim_onsite + M.Col[i]];
      }
      else {
         std::cout << "Error in DMRG_Make_Elem_Onsite_LLLR_RRRL at 1" << std::endl;
         std::exit(0);
      }
      
      if (inv < 0) {
         std::cout << "Error in DMRG_Make_Elem_Onsite_LLLR_RRRL at 2" << std::endl;
         std::cout << "Mat_Type=" << Mat_Type << std::endl;
         std::exit(0);
      }
      
      int check = A_Basis.Check_Basis[inv];
      
      if (check == -1) {
         A_Basis.Check_Basis[inv] = elem_num;
         A_Basis.Val[elem_num]    = coeff*M.Val[i];
         A_Basis.Inv[elem_num]    = inv;
         elem_num++;
      }
      else {
         A_Basis.Val[check] += coeff*M.Val[i];
      }
      
   }
   
   A_Basis.elem_num = elem_num;

}
