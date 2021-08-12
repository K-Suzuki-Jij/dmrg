//
//  Get_Ham_RRRL.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/12/12.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Ham_RRRL(Block_Operator &Enviro, DMRG_Basis_RRRL &Basis_RRRL, std::vector<int> &Ele_RR, int RR_site, Model_1D_AKLM &Model, CRS &Ham_RRRL) {

   int dim_RRRL   = (int)Basis_RRRL.RR.size();
   int dim_onsite = Model.dim_onsite;
   
   std::vector<long> Row_Elem_Num(dim_RRRL + 1, 0);
   std::vector<DMRG_A_Basis_Set> A_Basis;
   DMRG_Get_A_Basis_Set(A_Basis, dim_RRRL, Model.p_thread);
   
   for (int thread_num = 0; thread_num < Model.p_thread; thread_num++) {
      A_Basis[thread_num].base_RL        = dim_onsite;
      A_Basis[thread_num].RR_site        = RR_site;
      A_Basis[thread_num].zero_precision = Model.zero_precision;
   }
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_RRRL; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      Onsite_Basis.RR = Basis_RRRL.RR[i];
      Onsite_Basis.RL = Basis_RRRL.RL[i];

      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;

      Make_Elem_Ham_RRRL(Onsite_Basis, A_Basis[thread_num], Basis_RRRL.Inv, Enviro, Ele_RR, Model);
      
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
   for (int i = 0; i <= dim_RRRL; i++) {
      tot_elem_num += Row_Elem_Num[i];
   }
   
   //Do not use openmp here
   for (int i = 0; i < dim_RRRL; i++) {
      Row_Elem_Num[i + 1] += Row_Elem_Num[i];
   }
   
   Ham_RRRL.Col.resize(tot_elem_num);
   Ham_RRRL.Val.resize(tot_elem_num);
   Ham_RRRL.Row.resize(dim_RRRL + 1);

#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_RRRL; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      Onsite_Basis.RR = Basis_RRRL.RR[i];
      Onsite_Basis.RL = Basis_RRRL.RL[i];
      
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;

      Make_Elem_Ham_RRRL(Onsite_Basis, A_Basis[thread_num], Basis_RRRL.Inv, Enviro, Ele_RR, Model);
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int inv = A_Basis[thread_num].Inv[j];
         if (0 <= inv) {
            Ham_RRRL.Col[Row_Elem_Num[i]] = inv;
            Ham_RRRL.Val[Row_Elem_Num[i]] = A_Basis[thread_num].Val[j];
            Row_Elem_Num[i]++;
         }
         else {
            printf("Error in Get_Ham_LLLR at 2\n");
            exit(1);
         }
      }
      
      Ham_RRRL.Row[i + 1] = Row_Elem_Num[i];
      
   }
   
   if (Ham_RRRL.Row[dim_RRRL] != tot_elem_num) {
      std::cout << "Error in Get_Ham_LLLR at 3" << std::endl;
      exit(1);
   }
   
   Ham_RRRL.col_dim = dim_RRRL;
   Ham_RRRL.row_dim = dim_RRRL;
   
   Sort_Col_CRS(Ham_RRRL, Model.p_thread);

}
