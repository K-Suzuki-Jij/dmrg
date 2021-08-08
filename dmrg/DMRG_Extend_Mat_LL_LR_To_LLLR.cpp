#include <iostream>
#include "DMRG.hpp"

void DMRG_Extend_Mat_LL_LR_To_LLLR(const CRS &M_LL, const CRS &M_LR, CRS &Out_LLLR, const std::vector<int> &Ele_LL, std::string Sign_Flag, const DMRG_Basis_LLLR &Basis_LLLR) {
   
   int dim_LLLR   = (int)Basis_LLLR.LL.size();
   int dim_onsite = M_LR.row_dim;
   
   int c1 = (Sign_Flag == "LL_LR");
   int c2 = (Sign_Flag == "LR_LL");
   int c3 = (Sign_Flag == "No");
   
   int sign;
   
   Clear_Matrix(Out_LLLR);
   
   Out_LLLR.Row.push_back(0);
   
   for (int i = 0; i < dim_LLLR; i++) {
      int LL = Basis_LLLR.LL[i];
      int LR = Basis_LLLR.LR[i];
      
      if (c1) {
         if (Ele_LL[LL]%2 == 0) {
            sign = -1;
         }
         else {
            sign = 1;
         }
      }
      else if (c2) {
         if (Ele_LL[LL]%2 == 0) {
            sign = 1;
         }
         else {
            sign = -1;
         }
      }
      else if (c3) {
         sign = 1;
      }
      else {
         std::cout << "Error in DMRG_Extend_Mat_LL_LR_To_LLLR" << std::endl;
         exit(1);
      }
      
      for (long j = M_LL.Row[LL]; j < M_LL.Row[LL+1]; j++) {
         int    col_LL = M_LL.Col[j];
         double val_LL = M_LL.Val[j];
         for (long k = M_LR.Row[LR]; k < M_LR.Row[LR+1]; k++) {
            int    col_LR = M_LR.Col[k];
            double val_LR = M_LR.Val[k];
            int    inv    = Basis_LLLR.Inv[col_LL*dim_onsite + col_LR];
            if (inv >= 0) {
               Out_LLLR.Val.push_back(sign*val_LL*val_LR);
               Out_LLLR.Col.push_back(inv);
            }
         }
      }
      Out_LLLR.Row.push_back(Out_LLLR.Col.size());
   }
   
   Out_LLLR.row_dim = dim_LLLR;
   Out_LLLR.col_dim = dim_LLLR;
   
}
