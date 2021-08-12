//
//  Get_Ham_LRRL.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/11/16.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Ham_LRRL(DMRG_Basis_LRRL &Basis_LRRL, Model_1D_XXZ &Model, CRS &Ham_LRRL) {
   
   int dim_LRRL   = (int)Basis_LRRL.LR.size();
   int dim_onsite = Model.dim_onsite;
   
   std::vector<long> Row_Elem_Num(dim_LRRL + 1, 0);
   
   std::vector<DMRG_A_Basis_Set> A_Basis;
   DMRG_Get_A_Basis_Set(A_Basis, dim_LRRL, Model.p_thread);
   
   for (int thread_num = 0; thread_num < Model.p_thread; thread_num++) {
      A_Basis[thread_num].base_RL        = dim_onsite;
      A_Basis[thread_num].zero_precision = Model.zero_precision;
   }
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_LRRL; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      Onsite_Basis.LR = Basis_LRRL.LR[i];
      Onsite_Basis.RL = Basis_LRRL.RL[i];

      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;

      Make_Elem_Ham_LRRL(Onsite_Basis, A_Basis[thread_num], Basis_LRRL.Inv, Model);
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int inv = A_Basis[thread_num].Inv[j];
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
   for (int i = 0; i <= dim_LRRL; i++) {
      tot_elem_num += Row_Elem_Num[i];
   }
   
   //Do not use openmp here
   for (int i = 0; i < dim_LRRL; i++) {
      Row_Elem_Num[i + 1] += Row_Elem_Num[i];
   }
   
   Ham_LRRL.Col.resize(tot_elem_num);
   Ham_LRRL.Val.resize(tot_elem_num);
   Ham_LRRL.Row.resize(dim_LRRL + 1);
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_LRRL; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      Onsite_Basis.LR = Basis_LRRL.LR[i];
      Onsite_Basis.RL = Basis_LRRL.RL[i];
      
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;

      Make_Elem_Ham_LRRL(Onsite_Basis, A_Basis[thread_num], Basis_LRRL.Inv, Model);

      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int inv = A_Basis[thread_num].Inv[j];
         if (0 <= inv) {
            Ham_LRRL.Col[Row_Elem_Num[i]] = inv;
            Ham_LRRL.Val[Row_Elem_Num[i]] = A_Basis[thread_num].Val[j];
            Row_Elem_Num[i]++;
         }
         else {
            printf("Error in Get_Ham_LLLR at 2\n");
            exit(1);
         }
      }
      
      Ham_LRRL.Row[i + 1] = Row_Elem_Num[i];
      
   }
   
   if (Ham_LRRL.Row[dim_LRRL] != tot_elem_num) {
      std::cout << "Error in Get_Ham_LLLR at 3" << std::endl;
      exit(1);
   }
   
   Ham_LRRL.col_dim = dim_LRRL;
   Ham_LRRL.row_dim = dim_LRRL;
   
   Sort_Col_CRS(Ham_LRRL, Model.p_thread);
   
}

