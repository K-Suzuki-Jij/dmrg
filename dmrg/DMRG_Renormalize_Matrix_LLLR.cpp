#include <iostream>
#include "DMRG.hpp"

void DMRG_Renormalize_Matrix_LLLR(const CRS &M_LL, const CRS &M_LR, const CCS &Trans_Mat, CRS &M_Out_LL, const std::vector<int> &Ele_LL, std::string Sign_Flag, const DMRG_Basis_LLLR &Basis_LLLR) {
   
   if (Trans_Mat.col_dim > (int)Basis_LLLR.LL.size() || (int)Basis_LLLR.LL.size() != Trans_Mat.row_dim) {
      std::cout << "Error in DMRG_Renormalize_Matrix_LLLR" << std::endl;
      exit(1);
   }
   
   CRS M_LLLR;
   DMRG_Extend_Mat_LL_LR_To_LLLR(M_LL, M_LR, M_LLLR, Ele_LL, Sign_Flag, Basis_LLLR);
   
   Change_Matrix_Basis(M_LLLR, Trans_Mat, M_Out_LL, 1);
   
}
