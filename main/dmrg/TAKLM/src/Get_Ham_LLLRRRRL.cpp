//
//  Get_Ham_LLLRRRRL.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/05.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, DMRG_Time &Time, Model_1D_TAKLM &Model) {
   
   double start = omp_get_wtime();
   
   int LL_site           = Block.LL_site;
   int RR_site           = Block.RR_site;
   int dim_onsite        = (int)Basis_System.QN1_LL_LL_Stored[0].size();
   int dim_LL            = (int)Basis_System.QN1_LL_LL_Stored[LL_site].size();
   int dim_RR            = (int)Basis_Enviro.QN1_LL_LL_Stored[RR_site].size();
   int dim_LLLRRRRL      = (int)Basis.LLLRRRRL.LL.size();
   double zero_precision = Model.zero_precision;
   
   DMRG_Block_Hamiltonian Block_Ham;
   
   Block_Ham.Ele_LL.resize(dim_LL);
   Block_Ham.Ele_RR.resize(dim_RR);
   Block_Ham.Ele_On.resize(dim_onsite);
   
   for (int i = 0; i < dim_LL; i++) {
      Block_Ham.Ele_LL[i] = Basis_System.QN1_LL_LL_Stored[LL_site][i] + Basis_System.QN2_LL_LL_Stored[LL_site][i];
   }
   for (int i = 0; i < dim_RR; i++) {
      Block_Ham.Ele_RR[i] = Basis_Enviro.QN1_LL_LL_Stored[RR_site][i] + Basis_Enviro.QN2_LL_LL_Stored[RR_site][i];
   }
   for (int i = 0; i < dim_onsite; i++) {
      Block_Ham.Ele_On[i] = Basis_System.QN1_LL_LL_Stored[   0   ][i] + Basis_System.QN2_LL_LL_Stored[   0   ][i];
   }
   
   Get_Ham_Block("LLLR", System          , Basis, Block_Ham.Ele_LL, LL_site, int()  , Model, Block_Ham.LLLR);
   Get_Ham_Block("RRRL", Enviro          , Basis, Block_Ham.Ele_RR, int()  , RR_site, Model, Block_Ham.RRRL);
   Get_Ham_Block("LRRL", Block_Operator(), Basis, Block_Ham.Ele_On, int()  , RR_site, Model, Block_Ham.LRRL);
   
   std::vector<long> Row_Elem_Num(dim_LLLRRRRL + 1, 0);
   
   std::vector<DMRG_A_Basis_Set> A_Basis;
   DMRG_Get_A_Basis_Set(A_Basis, dim_LLLRRRRL, Model.p_thread);
   
   for (int thread_num = 0; thread_num < Model.p_thread; thread_num++) {
      A_Basis[thread_num].LL_site        = LL_site;
      A_Basis[thread_num].RR_site        = RR_site;
      A_Basis[thread_num].base_LRRRRL    = Basis.LLLRRRRL.base_LRRRRL;
      A_Basis[thread_num].base_RRRL      = Basis.LLLRRRRL.base_RRRL;
      A_Basis[thread_num].base_RL        = Basis.LLLRRRRL.base_RL;
      A_Basis[thread_num].zero_precision = zero_precision;
      A_Basis[thread_num].c_obc          = (Model.BC == "OBC");
   }
   
#pragma omp parallel for schedule(guided) num_threads (Model.p_thread)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      DMRG_Onsite_Basis Basis_Onsite;
      Basis_Onsite.row  = i;
      Basis_Onsite.LL   = Basis.LLLRRRRL.LL[i];
      Basis_Onsite.LR   = Basis.LLLRRRRL.LR[i];
      Basis_Onsite.RR   = Basis.LLLRRRRL.RR[i];
      Basis_Onsite.RL   = Basis.LLLRRRRL.RL[i];
      Basis_Onsite.LLLR = Basis.LLLR.Inv[Basis_Onsite.LL*dim_onsite + Basis_Onsite.LR];
      Basis_Onsite.RRRL = Basis.RRRL.Inv[Basis_Onsite.RR*dim_onsite + Basis_Onsite.RL];
      Basis_Onsite.LRRL = Basis.LRRL.Inv[Basis_Onsite.LR*dim_onsite + Basis_Onsite.RL];
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;
      DMRG_Make_Elem_Ham_LLLRRRRL(Basis_Onsite, A_Basis[thread_num], Basis, Block_Ham, "Yes");
      Row_Elem_Num[i + 1] += A_Basis[thread_num].elem_num;
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
   
   
#pragma omp parallel for schedule(guided) num_threads (Model.p_thread)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      DMRG_Onsite_Basis Basis_Onsite;
      Basis_Onsite.row  = i;
      Basis_Onsite.LL   = Basis.LLLRRRRL.LL[i];
      Basis_Onsite.LR   = Basis.LLLRRRRL.LR[i];
      Basis_Onsite.RR   = Basis.LLLRRRRL.RR[i];
      Basis_Onsite.RL   = Basis.LLLRRRRL.RL[i];
      Basis_Onsite.LLLR = Basis.LLLR.Inv[Basis_Onsite.LL*dim_onsite + Basis_Onsite.LR];
      Basis_Onsite.RRRL = Basis.RRRL.Inv[Basis_Onsite.RR*dim_onsite + Basis_Onsite.RL];
      Basis_Onsite.LRRL = Basis.LRRL.Inv[Basis_Onsite.LR*dim_onsite + Basis_Onsite.RL];
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;
      
      DMRG_Make_Elem_Ham_LLLRRRRL(Basis_Onsite, A_Basis[thread_num], Basis, Block_Ham, "Yes");
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int inv       = A_Basis[thread_num].Inv[j];
         Ham_LLLRRRRL.Col[Row_Elem_Num[i]] = inv;
         Ham_LLLRRRRL.Val[Row_Elem_Num[i]] = A_Basis[thread_num].Val[j];
         Row_Elem_Num[i]++;
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
   
   Time.make_ham = omp_get_wtime() - start;
   
}
