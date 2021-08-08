//
//  Created by Kohei Suzuki on 2021/01/07.
//

#include <iostream>
#include <cmath>
#include <omp.h>
#include "DMRG.hpp"

void DMRG_Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, const DMRG_Block_Hamiltonian &Block_Ham, const DMRG_Basis &Basis, const DMRG_Block_Information &Block, std::string Sign_Flag, int p_threads) {
   
   int LL_site      = Block.LL_site;
   int RR_site      = Block.RR_site;
   int dim_onsite   = Block.dim_onsite;
   int dim_RR       = Block.dim_RR;
   int dim_LLLRRRRL = (int)Basis.LLLRRRRL.LL.size();
   std::vector<long> Row_Elem_Num(dim_LLLRRRRL + 1, 0);
   
   std::vector<DMRG_A_Basis_Set> A_Basis;
   DMRG_Get_A_Basis_Set(A_Basis, dim_LLLRRRRL, p_threads);
   
   for (int thread_num = 0; thread_num < p_threads; thread_num++) {
      A_Basis[thread_num].LL_site        = LL_site;
      A_Basis[thread_num].RR_site        = RR_site;
      A_Basis[thread_num].base_LRRRRL    = Basis.LLLRRRRL.base_LRRRRL;
      A_Basis[thread_num].base_RRRL      = Basis.LLLRRRRL.base_RRRL;
      A_Basis[thread_num].base_RL        = Basis.LLLRRRRL.base_RL;
      A_Basis[thread_num].zero_precision = std::pow(10,-15);
   }
   
#pragma omp parallel for schedule(auto) num_threads (p_threads)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      DMRG_Onsite_Basis Basis_Onsite;
      Basis_Onsite.row  = i;
      Basis_Onsite.LL   = Basis.LLLRRRRL.LL[i];
      Basis_Onsite.LR   = Basis.LLLRRRRL.LR[i];
      Basis_Onsite.RR   = Basis.LLLRRRRL.RR[i];
      Basis_Onsite.RL   = Basis.LLLRRRRL.RL[i];
      Basis_Onsite.LLLR = Basis.LLLR.Inv[Basis_Onsite.LL*dim_onsite + Basis_Onsite.LR];
      Basis_Onsite.LLRR = Basis.LLRR.Inv[Basis_Onsite.LL*dim_RR     + Basis_Onsite.RR];
      Basis_Onsite.LLRL = Basis.LLRL.Inv[Basis_Onsite.LL*dim_onsite + Basis_Onsite.RL];
      Basis_Onsite.LRRR = Basis.LRRR.Inv[Basis_Onsite.LR*dim_RR     + Basis_Onsite.RR];
      Basis_Onsite.LRRL = Basis.LRRL.Inv[Basis_Onsite.LR*dim_onsite + Basis_Onsite.RL];
      Basis_Onsite.RRRL = Basis.RRRL.Inv[Basis_Onsite.RR*dim_onsite + Basis_Onsite.RL];
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;
      DMRG_Make_Elem_Ham_LLLRRRRL(Basis_Onsite, A_Basis[thread_num], Basis, Block_Ham, Sign_Flag);
      Row_Elem_Num[i + 1] += A_Basis[thread_num].elem_num;
   }
   
   long tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (p_threads)
   for (int i = 0; i <= dim_LLLRRRRL; i++) {
      tot_elem_num += Row_Elem_Num[i];
   }
   
   //Do not use openmp here
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      Row_Elem_Num[i + 1] += Row_Elem_Num[i];
   }
   
   Ham_LLLRRRRL.Col.resize(tot_elem_num);
   Ham_LLLRRRRL.Val.resize(tot_elem_num);
   Ham_LLLRRRRL.Row.resize(dim_LLLRRRRL + 1);
   
   
#pragma omp parallel for schedule(auto) num_threads (p_threads)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      DMRG_Onsite_Basis Basis_Onsite;
      Basis_Onsite.row  = i;
      Basis_Onsite.LL   = Basis.LLLRRRRL.LL[i];
      Basis_Onsite.LR   = Basis.LLLRRRRL.LR[i];
      Basis_Onsite.RR   = Basis.LLLRRRRL.RR[i];
      Basis_Onsite.RL   = Basis.LLLRRRRL.RL[i];
      Basis_Onsite.LLLR = Basis.LLLR.Inv[Basis_Onsite.LL*dim_onsite + Basis_Onsite.LR];
      Basis_Onsite.LLRR = Basis.LLRR.Inv[Basis_Onsite.LL*dim_RR     + Basis_Onsite.RR];
      Basis_Onsite.LLRL = Basis.LLRL.Inv[Basis_Onsite.LL*dim_onsite + Basis_Onsite.RL];
      Basis_Onsite.LRRR = Basis.LRRR.Inv[Basis_Onsite.LR*dim_RR     + Basis_Onsite.RR];
      Basis_Onsite.LRRL = Basis.LRRL.Inv[Basis_Onsite.LR*dim_onsite + Basis_Onsite.RL];
      Basis_Onsite.RRRL = Basis.RRRL.Inv[Basis_Onsite.RR*dim_onsite + Basis_Onsite.RL];
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;
      
      DMRG_Make_Elem_Ham_LLLRRRRL(Basis_Onsite, A_Basis[thread_num], Basis, Block_Ham, Sign_Flag);
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int inv = A_Basis[thread_num].Inv[j];
         Ham_LLLRRRRL.Col[Row_Elem_Num[i]] = inv;
         Ham_LLLRRRRL.Val[Row_Elem_Num[i]] = A_Basis[thread_num].Val[j];
         Row_Elem_Num[i]++;
      }
      
      Ham_LLLRRRRL.Row[i + 1] = Row_Elem_Num[i];
      
   }
   
   if (Ham_LLLRRRRL.Row[dim_LLLRRRRL] != tot_elem_num) {
      std::cout << "Error in DMRG_Get_Ham_LLLRRRRL at 1" << std::endl;
      std::exit(0);
   }
   
   Ham_LLLRRRRL.col_dim = dim_LLLRRRRL;
   Ham_LLLRRRRL.row_dim = dim_LLLRRRRL;
   
   Sort_Col_CRS(Ham_LLLRRRRL, p_threads);
   
}
