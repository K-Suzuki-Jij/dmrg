//
//  Get_Ham_Block.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/12/31.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Ham_Block(std::string Mat_Type, const Block_Operator &Block_Op, const DMRG_Basis &Basis, const std::vector<int> &Ele, int LL_site, int RR_site, const Model_1D_TAKLM &Model, CRS &Ham) {
   
   int dim_onsite = Model.dim_onsite;

   int c_LLLR = (Mat_Type == "LLLR");
   int c_LRRL = (Mat_Type == "LRRL");
   int c_RRRL = (Mat_Type == "RRRL");
   
   int dim_block = 0;
   const std::vector<int> *Inv;
        if (c_LLLR) {dim_block = (int)Basis.LLLR.LL.size(); Inv = &Basis.LLLR.Inv;}
   else if (c_LRRL) {dim_block = (int)Basis.LRRL.LR.size(); Inv = &Basis.LRRL.Inv;}
   else if (c_RRRL) {dim_block = (int)Basis.RRRL.RR.size(); Inv = &Basis.RRRL.Inv;}
   else             {std::cout << "Error in  Get_Ham_Block at 1" << std::endl; std::exit(0);}
   
   std::vector<long> Row_Elem_Num(dim_block + 1, 0);
   std::vector<DMRG_A_Basis_Set> A_Basis;
   
   DMRG_Get_A_Basis_Set(A_Basis, dim_block, Model.p_thread);
   
   for (int thread_num = 0; thread_num < Model.p_thread; thread_num++) {
      A_Basis[thread_num].base_RL        = dim_onsite;
      A_Basis[thread_num].LL_site        = LL_site;
      A_Basis[thread_num].RR_site        = RR_site;
      A_Basis[thread_num].zero_precision = Model.zero_precision;
   }
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_block; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      
      if (c_LLLR) {
         Onsite_Basis.LL = Basis.LLLR.LL[i];
         Onsite_Basis.LR = Basis.LLLR.LR[i];
      }
      else if (c_LRRL) {
         Onsite_Basis.LR = Basis.LRRL.LR[i];
         Onsite_Basis.RL = Basis.LRRL.RL[i];
      }
      else if (c_RRRL) {
         Onsite_Basis.RR = Basis.RRRL.RR[i];
         Onsite_Basis.RL = Basis.RRRL.RL[i];
      }
      else {
         std::cout << "Error in  Get_Ham_Block at 2" << std::endl;
         std::exit(0);
      }

      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;

      Make_Elem_Ham(Mat_Type, Onsite_Basis, A_Basis[thread_num], *Inv, Block_Op, Ele, Model);
      
      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         if (0 <= A_Basis[thread_num].Inv[j]) {
            Row_Elem_Num[i + 1]++;
         }
         else {
            std::cout << "Error in  Get_Ham_Block at 3" << std::endl;
            std::exit(0);
         }
      }
   }
   
   long tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (Model.p_thread)
   for (int i = 0; i <= dim_block; i++) {
      tot_elem_num += Row_Elem_Num[i];
   }
   
   //Do not use openmp here
   for (int i = 0; i < dim_block; i++) {
      Row_Elem_Num[i + 1] += Row_Elem_Num[i];
   }
   
   Ham.Col.resize(tot_elem_num);
   Ham.Val.resize(tot_elem_num);
   Ham.Row.resize(dim_block + 1);
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_block; i++) {
      DMRG_Onsite_Basis Onsite_Basis;
      
      if (c_LLLR) {
         Onsite_Basis.LL = Basis.LLLR.LL[i];
         Onsite_Basis.LR = Basis.LLLR.LR[i];
      }
      else if (c_LRRL) {
         Onsite_Basis.LR = Basis.LRRL.LR[i];
         Onsite_Basis.RL = Basis.LRRL.RL[i];
      }
      else if (c_RRRL) {
         Onsite_Basis.RR = Basis.RRRL.RR[i];
         Onsite_Basis.RL = Basis.RRRL.RL[i];
      }
      else {
         std::cout << "Error in  Get_Ham_Block at 2" << std::endl;
         std::exit(0);
      }
      
      int thread_num = omp_get_thread_num();
      A_Basis[thread_num].elem_num = 0;

      Make_Elem_Ham(Mat_Type, Onsite_Basis, A_Basis[thread_num], *Inv, Block_Op, Ele, Model);

      for (int j = 0; j < A_Basis[thread_num].elem_num; j++) {
         int inv = A_Basis[thread_num].Inv[j];
         if (0 <= inv) {
            Ham.Col[Row_Elem_Num[i]] = inv;
            Ham.Val[Row_Elem_Num[i]] = A_Basis[thread_num].Val[j];
            Row_Elem_Num[i]++;
         }
         else {
            std::cout << "Error in  Get_Ham_Block at 3" << std::endl;
            std::exit(0);
         }
      }
      
      Ham.Row[i + 1] = Row_Elem_Num[i];
   }
   
   if (Ham.Row[dim_block] != tot_elem_num) {
      std::cout << "Error in Get_Ham_LLLR at 4" << std::endl;
      exit(1);
   }
   
   Ham.col_dim = dim_block;
   Ham.row_dim = dim_block;
   
   Sort_Col_CRS(Ham, Model.p_thread);
   
   
}
