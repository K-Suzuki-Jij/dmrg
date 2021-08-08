//
//  Created by Kohei Suzuki on 2020/12/28.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include <iostream>
#include "DMRG.hpp"

void DMRG_Make_Elem_Onsite_LLLRRRRL(const std::string Mat_Type, int basis_change, const DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, const std::vector<short> &Basis_1, const std::vector<short> &Basis_2, const std::vector<int> &Inv_LLLRRRRL, const CRS &M, double coeff, const std::vector<int> &Ele_RR, const std::vector<int> &Ele_On, const std::string Sign_Flag) {
   
   if (std::abs(coeff) < A_Basis.zero_precision) {
      return;
   }
   
   int LL          = Basis.LL;
   int LR          = Basis.LR;
   int RR          = Basis.RR;
   int RL          = Basis.RL;
   int row         = Basis.row;
   int elem_num    = A_Basis.elem_num;
   int base_LRRRRL = A_Basis.base_LRRRRL;
   int base_RRRL   = A_Basis.base_RRRL;
   int base_RL     = A_Basis.base_RL;
   int sign_col    = 1;
   
   int c_LLLR = (Mat_Type == "LLLR");//  No Sign
   int c_LLRR = (Mat_Type == "LLRR");//Need Sign
   int c_LLRL = (Mat_Type == "LLRL");//Need Sign
   int c_LRRR = (Mat_Type == "LRRR");//  No Sign
   int c_LRRL = (Mat_Type == "LRRL");//Need Sign
   int c_RRRL = (Mat_Type == "RRRL");//  No Sign
   int c_sign = (Sign_Flag == "Yes");
   
   if (c_sign) {
      if ((c_LLRR || c_LRRL) && Ele_On[LR]%2 == 1 && Ele_RR[RR]%2 == 1) {
         sign_col = -1;
      }
      if (c_LLRL && (Ele_On[LR] + Ele_RR[RR])%2 == 1 && Ele_On[RL]%2 == 1) {
         sign_col = -1;
      }
   }

   
   for (long i = M.Row[basis_change]; i < M.Row[basis_change + 1]; i++) {
      int col_1 = Basis_1[M.Col[i]];
      int col_2 = Basis_2[M.Col[i]];
      int inv = -1;
      if (c_LLLR) {
         inv = Inv_LLLRRRRL[(long)col_1*base_LRRRRL + col_2*base_RRRL + RR*base_RL + RL];
      }
      else if (c_LLRR) {
         inv = Inv_LLLRRRRL[(long)col_1*base_LRRRRL + LR*base_RRRL + col_2*base_RL + RL];
      }
      else if (c_LLRL) {
         inv = Inv_LLLRRRRL[(long)col_1*base_LRRRRL + LR*base_RRRL + RR*base_RL + col_2];
      }
      else if (c_LRRR) {
         inv = Inv_LLLRRRRL[(long)LL*base_LRRRRL + col_1*base_RRRL + col_2*base_RL + RL];
      }
      else if (c_LRRL) {
         inv = Inv_LLLRRRRL[(long)LL*base_LRRRRL + col_1*base_RRRL + RR*base_RL + col_2];
      }
      else if (c_RRRL) {
         inv = Inv_LLLRRRRL[(long)LL*base_LRRRRL + LR*base_RRRL + col_1*base_RL + col_2];
      }
      else {
         std::cout << "Error in DMRG_Make_Elem_Onsite_LLLRRRRL at 1" << std::endl;
         std::exit(0);
      }
      
      if (inv < 0) {
         std::cout << "Error in DMRG_Make_Elem_Onsite_LLLRRRRL at 2" << std::endl;
         std::cout << Mat_Type << ", row=" << basis_change << ", col=" << M.Col[i] << ", inv=" << inv << ", val=" << M.Val[i] << std::endl;
         std::exit(0);
      }
      
      int sign_row = 1;

      if (0 <= inv && inv <= row) {
         if (c_sign) {
            if (c_LLRR && Ele_On[LR]%2 == 1 && Ele_RR[col_2]%2 == 1) {
               sign_row = -1;
            }
            if (c_LRRL && Ele_On[col_1]%2 == 1 && Ele_RR[RR]%2 == 1) {
               sign_row = -1;
            }
            if (c_LLRL && (Ele_On[LR] + Ele_RR[RR])%2 == 1 && Ele_On[col_2]%2 == 1) {
               sign_row = -1;
            }
         }

         int check = A_Basis.Check_Basis[inv];
         if (check == -1) {
            A_Basis.Check_Basis[inv] = elem_num;
            A_Basis.Val[elem_num]    = coeff*sign_row*sign_col*M.Val[i];
            A_Basis.Inv[elem_num]    = inv;
            elem_num++;
         }
         else {
            A_Basis.Val[check] += coeff*sign_row*sign_col*M.Val[i];
         }
      }
   }
   
   A_Basis.elem_num = elem_num;
   
}
