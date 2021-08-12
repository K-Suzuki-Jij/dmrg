//
//  Get_Ham_LLLR.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Ham_LLLR(Block_Operator &System, DMRG_LLLR_Basis &Basis_LLLR, DMRG_Block_Composition &Block_Compo, std::vector<int> &Ele_LL, Model_1D_AKLM_TVF &Model, CRS &Ham_LLLR) {
   
   int dim_LLLR   = Block_Compo.dim_LLLR;
   int dim_onsite = Block_Compo.dim_onsite;
   
   std::vector<long> Row_Elem_Num(dim_LLLR + 1, 0);
   
   std::vector<DMRG_A_Basis_Set> A_Basis;
   DMRG_Get_A_Basis_Set(A_Basis, dim_LLLR, Model.p_thread);
   
   for (int thread_num = 0; thread_num < Model.p_thread; thread_num++) {
      A_Basis[thread_num].base_RL        = dim_onsite;
      A_Basis[thread_num].LL_site        = Block_Compo.LL_site;
      A_Basis[thread_num].zero_precision = Model.zero_precision;
   }
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_LLLR; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      Onsite_Basis.LL = Basis_LLLR.LL[i];
      Onsite_Basis.LR = Basis_LLLR.LR[i];

      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;

      Make_Elem_Ham_LLLR(Onsite_Basis, A_Basis[thread_num], Basis_LLLR.Inv, System, Ele_LL, Model);
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int LL  = A_Basis[thread_num].LL[j];
         int LR  = A_Basis[thread_num].LR[j];
         int inv = Basis_LLLR.Inv[LL*dim_onsite + LR];
         if (0 <= inv) {
            Row_Elem_Num[i + 1]++;
         }
         else {
            printf("Error in Get_Ham_LLLR at 1\n");
            exit(1);
         }
      }
   }
   
   long tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (Model.p_thread)
   for (int i = 0; i <= dim_LLLR; i++) {
      tot_elem_num += Row_Elem_Num[i];
   }
   
   //Do not use openmp here
   for (int i = 0; i < dim_LLLR; i++) {
      Row_Elem_Num[i + 1] += Row_Elem_Num[i];
   }
   
   Ham_LLLR.Col.resize(tot_elem_num);
   Ham_LLLR.Val.resize(tot_elem_num);
   Ham_LLLR.Row.resize(dim_LLLR + 1);
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_LLLR; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      Onsite_Basis.LL = Basis_LLLR.LL[i];
      Onsite_Basis.LR = Basis_LLLR.LR[i];
      
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;

      Make_Elem_Ham_LLLR(Onsite_Basis, A_Basis[thread_num], Basis_LLLR.Inv, System, Ele_LL, Model);
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int LL  = A_Basis[thread_num].LL[j];
         int LR  = A_Basis[thread_num].LR[j];
         int inv = Basis_LLLR.Inv[LL*dim_onsite + LR];
         if (0 <= inv) {
            Ham_LLLR.Col[Row_Elem_Num[i]] = inv;
            Ham_LLLR.Val[Row_Elem_Num[i]] = A_Basis[thread_num].Val[j];
            Row_Elem_Num[i]++;
         }
         else {
            printf("Error in Get_Ham_LLLR at 2\n");
            exit(1);
         }
      }
      
      Ham_LLLR.Row[i + 1] = Row_Elem_Num[i];
      
   }
   
   if (Ham_LLLR.Row[dim_LLLR] != tot_elem_num) {
      std::cout << "Error in Get_Ham_LLLR at 3" << std::endl;
      exit(1);
   }
   
   Ham_LLLR.col_dim = dim_LLLR;
   Ham_LLLR.row_dim = dim_LLLR;
   
   Sort_Col_CRS(Ham_LLLR, Model.p_thread);
   
}
