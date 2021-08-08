#include <iostream>
#include "DMRG.hpp"

void DMRG_Clear_Check_Basis(DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv) {
   
   for (int i = 0; i < A_Basis.elem_num; i++) {
      int inv = A_Basis.Inv[i];
      A_Basis.Check_Basis[inv] = -1;
   }
   
}
