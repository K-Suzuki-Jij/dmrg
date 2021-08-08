#include <iostream>
#include <algorithm>
#include "DMRG.hpp"

//Quick sort
void DMRG_Sort_Bases(DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, std::vector<std::vector<int>> &QN_LLLR_LLLRRRRL, int left, int right) {

   if (right - left <= 1) {
      return;
   }
   
   int   num_of_qn   = (int)QN_LLLR_LLLRRRRL.size();
   int   pivot_index = (left + right)/2;
   short pivot_LL    = Basis_LLLRRRRL.LL[pivot_index];
   short pivot_LR    = Basis_LLLRRRRL.LR[pivot_index];
   short pivot_RR    = Basis_LLLRRRRL.RR[pivot_index];
   short pivot_RL    = Basis_LLLRRRRL.RL[pivot_index];
   
   int Pivot_QN[num_of_qn];
   for (int i = 0; i < num_of_qn; i++) {
      Pivot_QN[i] = QN_LLLR_LLLRRRRL[i][pivot_index];
   }
   
   std::swap(Basis_LLLRRRRL.LL[pivot_index], Basis_LLLRRRRL.LL[right - 1]);
   std::swap(Basis_LLLRRRRL.LR[pivot_index], Basis_LLLRRRRL.LR[right - 1]);
   std::swap(Basis_LLLRRRRL.RR[pivot_index], Basis_LLLRRRRL.RR[right - 1]);
   std::swap(Basis_LLLRRRRL.RL[pivot_index], Basis_LLLRRRRL.RL[right - 1]);
   
   for (int i = 0; i < num_of_qn; i++) {
      std::swap(QN_LLLR_LLLRRRRL[i][pivot_index], QN_LLLR_LLLRRRRL[i][right - 1]);
   }
   
   int i = left;
   
   for (int j = left; j < right - 1; j++) {
      int e_LL  = (Basis_LLLRRRRL.LL[j] == pivot_LL);
      int e_LR  = (Basis_LLLRRRRL.LR[j] == pivot_LR);
      int e_RR  = (Basis_LLLRRRRL.RR[j] == pivot_RR);
      int E_QN[num_of_qn];
      for (int k = 0; k < num_of_qn; k++) {
         E_QN[k] = (QN_LLLR_LLLRRRRL[k][j] == Pivot_QN[k]);
      }
      
      int c_qn     = 0;
      int c_qn_tot = 1;
      for (int k = 0; k < num_of_qn; k++) {
         c_qn_tot = (c_qn_tot && E_QN[k]);
         int c_qn_temp  = (QN_LLLR_LLLRRRRL[k][j] < Pivot_QN[k]);
         for (int l = 0; l < k; l++) {
            c_qn_temp = c_qn_temp && E_QN[l];
         }
         c_qn = c_qn || c_qn_temp;
      }
      
      int c_LL  = c_qn_tot && (Basis_LLLRRRRL.LL[j] < pivot_LL);
      int c_LR  = c_qn_tot && e_LL && (Basis_LLLRRRRL.LR[j] < pivot_LR);
      int c_RR  = c_qn_tot && e_LL && e_LR && (Basis_LLLRRRRL.RR[j] < pivot_RR);
      int c_RL  = c_qn_tot && e_LL && e_LR && e_RR && (Basis_LLLRRRRL.RL[j] < pivot_RL);
      
      if (c_qn || c_LL || c_LR || c_RR || c_RL) {
         std::swap(Basis_LLLRRRRL.LL[i], Basis_LLLRRRRL.LL[j]);
         std::swap(Basis_LLLRRRRL.LR[i], Basis_LLLRRRRL.LR[j]);
         std::swap(Basis_LLLRRRRL.RR[i], Basis_LLLRRRRL.RR[j]);
         std::swap(Basis_LLLRRRRL.RL[i], Basis_LLLRRRRL.RL[j]);
         for (int k = 0; k < num_of_qn; k++) {
            std::swap(QN_LLLR_LLLRRRRL[k][i], QN_LLLR_LLLRRRRL[k][j]);
         }
         i++;
      }
   }
   
   std::swap(Basis_LLLRRRRL.LL[i], Basis_LLLRRRRL.LL[right - 1]);
   std::swap(Basis_LLLRRRRL.LR[i], Basis_LLLRRRRL.LR[right - 1]);
   std::swap(Basis_LLLRRRRL.RR[i], Basis_LLLRRRRL.RR[right - 1]);
   std::swap(Basis_LLLRRRRL.RL[i], Basis_LLLRRRRL.RL[right - 1]);
   for (int k = 0; k < num_of_qn; k++) {
      std::swap(QN_LLLR_LLLRRRRL[k][i], QN_LLLR_LLLRRRRL[k][right - 1]);
   }
   
   DMRG_Sort_Bases(Basis_LLLRRRRL, QN_LLLR_LLLRRRRL, left, i    );
   DMRG_Sort_Bases(Basis_LLLRRRRL, QN_LLLR_LLLRRRRL, i+1 , right);

}
