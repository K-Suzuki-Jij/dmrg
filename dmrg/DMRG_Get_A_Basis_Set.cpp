#include "DMRG.hpp"

void DMRG_Get_A_Basis_Set(std::vector<DMRG_A_Basis_Set> &A_Basis, int num, int p_threads) {
   
   A_Basis.resize(p_threads);
   
   for (int thread_num = 0; thread_num < p_threads; thread_num++) {
      A_Basis[thread_num].Inv.resize(num);
      A_Basis[thread_num].Val.resize(num);
      A_Basis[thread_num].Check_Basis.resize(num);
   }
   
#pragma omp parallel for num_threads (p_threads)
   for (int thread_num = 0; thread_num < p_threads; thread_num++) {
      for (int j = 0; j < num; j++) {
         A_Basis[thread_num].Inv[j] = -1;
         A_Basis[thread_num].Val[j] = 0.0;
         A_Basis[thread_num].Check_Basis[j] = -1;
      }
   }
   
   
}
