#include <iostream>
#include "DMRG.hpp"

void DMRG_Make_Elem_Zero_LLLRRRRL(const DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv_LLLRRRRL) {

   long LLLRRRRL = (long)Basis.LL*A_Basis.base_LRRRRL + Basis.LR*A_Basis.base_RRRL + Basis.RR*A_Basis.base_RL + Basis.RL;
   int  inv      = Inv_LLLRRRRL[LLLRRRRL];
   if (inv >= 0) {
      int  check = A_Basis.Check_Basis[inv];
      if (check == -1) {
         A_Basis.Check_Basis[inv] = A_Basis.elem_num;
         A_Basis.Val[A_Basis.elem_num] = 0.0;
         A_Basis.Inv[A_Basis.elem_num] = inv;
         A_Basis.elem_num++;
      }
   }
   else {
      std::cout << "Error in DMRG_Make_Elem_zero_LLLRRRRL" << std::endl;
      exit(1);
   }
   
   
}
