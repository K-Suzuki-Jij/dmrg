//
//  Diagonalize_Ham_LLLRRRRL.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Diagonalize_Ham_LLLRRRRL(Block_Operator &System, Block_Operator &Enviro, DMRG_Basic_Information &Info, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, Model_1D_AKLM_TVF &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param, DMRG_Time &Time) {
   
   double start = omp_get_wtime();
   
   Diag_Param.Mat_Type = Model.Mat_Type;
   
   CRS Ham_LLLRRRRL;
   Get_Ham_LLLRRRRL(Ham_LLLRRRRL, Info.Basis_LLLRRRRL, System, Enviro, Basis_System, Basis_Enviro, Info.Block_Compo, Model);
   Time.make_ham = omp_get_wtime() - start;
   
   //Release Inv_LLLRRRRL
   std::vector<int>().swap(Info.Basis_LLLRRRRL.Inv);
   
   if (Dmrg_Param.Initial_Guess == "Yes" && Info.GS.Vector_Guess.size() == Ham_LLLRRRRL.row_dim) {
      Info.GS.Vector = Info.GS.Vector_Guess;
      Diag_Param.Lanczos_Initial_Guess = "Yes";
      Diag_Param.Calc_Vec = "No";
   }
   else {
      Diag_Param.Lanczos_Initial_Guess = "No";
      Diag_Param.Calc_Vec = "Yes";
   }
   
   start = omp_get_wtime();
   Lanczos(Ham_LLLRRRRL, Info.GS.Vector, Info.GS.val, Diag_Param, Model.p_thread);
   Time.diag = omp_get_wtime() - start;
   
   start = omp_get_wtime();
   Inverse_Iteration(Ham_LLLRRRRL, Info.GS.Vector, Info.GS.val, Diag_Param, Model.p_thread);
   Info.GS.error = Residual_Error_Eigenpair(Ham_LLLRRRRL, Info.GS.Vector, Info.GS.val, Diag_Param.Mat_Type, Model.p_thread);
   Time.inv_iter = omp_get_wtime() - start;
   
   //Reallocate Inv_LLLRRRRL
   DMRG_Get_Inv_LLLRRRRL(Info.Basis_LLLRRRRL.Inv, Info.Basis_LLLRRRRL, Info.Block_Compo, Model.p_thread);
   
}
