//
//  Diagonalize_Ham_LLLRRRRL.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/07/09.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Diagonalize_Ham_LLLRRRRL(DMRG_Ground_State &GS, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_AKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param, DMRG_Time &Time) {
   
   double start = omp_get_wtime();
   
   Diag_Param.Mat_Type = "Sym";
   
   CRS Ham_LLLRRRRL;
   Get_Ham_LLLRRRRL(Ham_LLLRRRRL, System, Enviro, Basis, Basis_System, Basis_Enviro, Block, Model);

   Time.make_ham = omp_get_wtime() - start;
   
   //Release Inv_LLLRRRRL
   std::vector<int>().swap(Basis.LLLRRRRL.Inv);
   
   if (Dmrg_Param.Initial_Guess == "Yes" && (int)GS.Vector_Guess.size() == Ham_LLLRRRRL.row_dim) {
      GS.Vector = GS.Vector_Guess;
      Diag_Param.Lanczos_Initial_Guess = "Yes";
      Diag_Param.Calc_Vec = "No";
   }
   else {
      Diag_Param.Lanczos_Initial_Guess = "No";
      Diag_Param.Calc_Vec = "Yes";
   }

   start = omp_get_wtime();
   Lanczos(Ham_LLLRRRRL, GS.Vector, GS.val, Diag_Param, Model.p_thread);
   Time.diag = omp_get_wtime() - start;
   
   start = omp_get_wtime();
   Inverse_Iteration(Ham_LLLRRRRL, GS.Vector, GS.val, Diag_Param, Model.p_thread);
   GS.error = Residual_Error_Eigenpair(Ham_LLLRRRRL, GS.Vector, GS.val, Diag_Param.Mat_Type, Model.p_thread);
   Time.inv_iter = omp_get_wtime() - start;
   
   //Reallocate Inv_LLLRRRRL
   int LL_site      = Block.LL_site;
   int RR_site      = Block.RR_site;
   int dim_LL       = (int)Basis_System.QN1_LL_LL_Stored[LL_site].size();
   int dim_RR       = (int)Basis_Enviro.QN1_LL_LL_Stored[RR_site].size();
   int dim_onsite   = (int)Basis_System.QN1_LL_LL_Stored[   0   ].size();
   int dim_LLLRRRRL = (int)Basis.LLLRRRRL.LL.size();
   int base_LRRRRL  = dim_onsite*dim_RR*dim_onsite;
   int base_RRRL    = dim_RR*dim_onsite;
   long whole_dim   = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   
   Basis.LLLRRRRL.Inv.resize(whole_dim);
#pragma omp parallel for num_threads (Model.p_thread)
   for (long i = 0; i < whole_dim; i++) {
      Basis.LLLRRRRL.Inv[i] = -1;
   }
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      int  LL       = Basis.LLLRRRRL.LL[i];
      int  LR       = Basis.LLLRRRRL.LR[i];
      int  RR       = Basis.LLLRRRRL.RR[i];
      int  RL       = Basis.LLLRRRRL.RL[i];
      long LLLRRRRL = (long)LL*base_LRRRRL + LR*base_RRRL + RR*dim_onsite + RL;
      Basis.LLLRRRRL.Inv[LLLRRRRL] = i;
   }
   
}
