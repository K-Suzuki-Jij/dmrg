//
//  Created by Kohei Suzuki on 2020/12/31.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include <omp.h>
#include "SML.hpp"
#include "DMRG.hpp"

void DMRG_Diagonalize_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, DMRG_Ground_State &GS, DMRG_Basis &Basis, DMRG_Block_Information &Block, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param, DMRG_Time &Time, int p_threads) {


   Diag_Param.Mat_Type = "Sym";
   Diag_Param.Calc_Vec = "Yes";

   //Release Inv_LLLRRRRL
   std::vector<int>().swap(Basis.LLLRRRRL.Inv);
   
   if (Dmrg_Param.Initial_Guess == "Yes" && (int)GS.Vector_Guess.size() == Ham_LLLRRRRL.row_dim) {
      GS.Vector = GS.Vector_Guess;
      Diag_Param.Lanczos_Initial_Guess = "Yes";
   }
   else {
      Diag_Param.Lanczos_Initial_Guess = "No";
   }

   double start = omp_get_wtime();
   Lanczos(Ham_LLLRRRRL, GS.Vector, GS.val, Diag_Param, p_threads);
   Time.diag = omp_get_wtime() - start;
   
   start = omp_get_wtime();
   Inverse_Iteration(Ham_LLLRRRRL, GS.Vector, GS.val, Diag_Param, p_threads);
   GS.error = Residual_Error_Eigenpair(Ham_LLLRRRRL, GS.Vector, GS.val, Diag_Param.Mat_Type, p_threads);
   Time.inv_iter = omp_get_wtime() - start;
   
   DMRG_Get_Inv_LLLRRRRL(Basis.LLLRRRRL, Block, p_threads);
   
}
