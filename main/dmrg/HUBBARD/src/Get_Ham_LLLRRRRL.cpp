//
//  Get_Ham_LLLRRRRL.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, DMRG_LLLRRRRL_Basis &Basis_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Composition &Block_Compo, Model_1D_HUBBARD &Model) {
   
   int LL_site           = Block_Compo.LL_site;
   int RR_site           = Block_Compo.RR_site;
   int dim_LLLRRRRL      = Block_Compo.dim_LLLRRRRL;
   int dim_onsite        = Model.dim_onsite;
   int Non_Sym_Condition = (Model.Mat_Type == "Non_Sym");
   int Sym_Condition     = (Model.Mat_Type == "Sym");
   int base_LRRRRL       = Basis_LLLRRRRL.base_LRRRRL;
   int base_RRRL         = Basis_LLLRRRRL.base_RRRL;
   double zero_precision = Model.zero_precision;
   
   std::vector<long> Row_Elem_Num(dim_LLLRRRRL + 1, 0);
   
   std::vector<DMRG_A_Basis_Set> A_Basis;
   DMRG_Get_A_Basis_Set(A_Basis, dim_LLLRRRRL, Model.p_thread);
   
   for (int thread_num = 0; thread_num < Model.p_thread; thread_num++) {
      A_Basis[thread_num].LL_site        = Block_Compo.LL_site;
      A_Basis[thread_num].RR_site        = Block_Compo.RR_site;
      A_Basis[thread_num].base_LRRRRL    = Basis_LLLRRRRL.base_LRRRRL;
      A_Basis[thread_num].base_RRRL      = Basis_LLLRRRRL.base_RRRL;
      A_Basis[thread_num].base_RL        = Basis_LLLRRRRL.base_RL;
      A_Basis[thread_num].zero_precision = zero_precision;
      A_Basis[thread_num].c_obc          = (Model.BC == "OBC");
   }
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      Onsite_Basis.LL = Basis_LLLRRRRL.LL[i];
      Onsite_Basis.LR = Basis_LLLRRRRL.LR[i];
      Onsite_Basis.RR = Basis_LLLRRRRL.RR[i];
      Onsite_Basis.RL = Basis_LLLRRRRL.RL[i];
      
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;
      
      Make_Elem_Ham_LLLRRRRL(Onsite_Basis,
                             A_Basis[thread_num],
                             Basis_LLLRRRRL.Inv,
                             System,
                             Enviro,
                             Basis_System.QN1_LL_LL_Stored[LL_site],
                             Basis_System.QN1_LL_LL_Stored[0],
                             Basis_Enviro.QN1_LL_LL_Stored[RR_site],
                             Model
                             );
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int LL        = A_Basis[thread_num].LL[j];
         int LR        = A_Basis[thread_num].LR[j];
         int RR        = A_Basis[thread_num].RR[j];
         int RL        = A_Basis[thread_num].RL[j];
         long LLLRRRRL = (long)LL*base_LRRRRL + LR*base_RRRL + RR*dim_onsite + RL;
         int inv       = Basis_LLLRRRRL.Inv[LLLRRRRL];
         int c_zero    = (std::fabs(A_Basis[thread_num].Val[j]) > zero_precision);
         int c1        = (0 <= inv) && (c_zero) && Non_Sym_Condition;
         int c2        = (0 <= inv) && (inv <= i) && (c_zero) && Sym_Condition;
         if (c1 || c2 || inv == i) {
            Row_Elem_Num[i + 1]++;
         }
         else if (inv < 0) {
            printf("Error in Get_Ham_LLLRRRRL at 1\n");
            exit(1);
         }
      }
   }
   
   long tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (Model.p_thread)
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
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      Onsite_Basis.LL = Basis_LLLRRRRL.LL[i];
      Onsite_Basis.LR = Basis_LLLRRRRL.LR[i];
      Onsite_Basis.RR = Basis_LLLRRRRL.RR[i];
      Onsite_Basis.RL = Basis_LLLRRRRL.RL[i];
      
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;
      
      Make_Elem_Ham_LLLRRRRL(Onsite_Basis,
                             A_Basis[thread_num],
                             Basis_LLLRRRRL.Inv,
                             System,
                             Enviro,
                             Basis_System.QN1_LL_LL_Stored[LL_site],
                             Basis_System.QN1_LL_LL_Stored[0],
                             Basis_Enviro.QN1_LL_LL_Stored[RR_site],
                             Model
                             );
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int LL        = A_Basis[thread_num].LL[j];
         int LR        = A_Basis[thread_num].LR[j];
         int RR        = A_Basis[thread_num].RR[j];
         int RL        = A_Basis[thread_num].RL[j];
         long LLLRRRRL = (long)LL*base_LRRRRL + LR*base_RRRL + RR*dim_onsite + RL;
         int inv       = Basis_LLLRRRRL.Inv[LLLRRRRL];
         int c_zero    = (std::fabs(A_Basis[thread_num].Val[j]) > zero_precision);
         int c1        = (0 <= inv) && (c_zero) && Non_Sym_Condition;
         int c2        = (0 <= inv) && (inv <= i) && (c_zero) && Sym_Condition;
         if (c1 || c2 || inv == i) {
            Ham_LLLRRRRL.Col[Row_Elem_Num[i]] = inv;
            Ham_LLLRRRRL.Val[Row_Elem_Num[i]] = A_Basis[thread_num].Val[j];
            Row_Elem_Num[i]++;
         }
         else if (inv < 0) {
            printf("Error in Get_Ham_LLLRRRRL at 2\n");
            exit(1);
         }
      }
      
      Ham_LLLRRRRL.Row[i + 1] = Row_Elem_Num[i];
      
   }
   
   if (Ham_LLLRRRRL.Row[dim_LLLRRRRL] != tot_elem_num) {
      printf("Error in Get_Ham_LLLRRRRL at 3\n");
      exit(1);
   }
   
   Ham_LLLRRRRL.col_dim = dim_LLLRRRRL;
   Ham_LLLRRRRL.row_dim = dim_LLLRRRRL;
   
   Sort_Col_CRS(Ham_LLLRRRRL, Model.p_thread);
   
}
