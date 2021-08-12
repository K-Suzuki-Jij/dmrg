//
//  Renormalize_System.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/06/19.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Renormalize_System(Block_Operator &System, DMRG_Basis_LLLR &Basis_LLLR, DMRG_Basis_Stored &Basis_System, int LL_site, DMRG_T_Mat &T_Mat, Model_1D_XXZ &Model) {
   
   CRS Ham_LLLR;
   Get_Ham_LLLR(System, Basis_LLLR, LL_site, Model, Ham_LLLR);
   
   Change_Matrix_Basis(Ham_LLLR, T_Mat.Trans_Mat, System.Ham[LL_site + 1]);
   
   std::vector<int> Dummy;
#pragma omp parallel sections num_threads (Model.p_thread)
   {
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.Sz_On, T_Mat.Trans_Mat, System.Sz_RE[LL_site + 1], Dummy, "No", Basis_LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.Sp_On, T_Mat.Trans_Mat, System.Sp_RE[LL_site + 1], Dummy, "No", Basis_LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.Sm_On, T_Mat.Trans_Mat, System.Sm_RE[LL_site + 1], Dummy, "No", Basis_LLLR);
   }
   
   int dim_renorm = T_Mat.Trans_Mat.col_dim;
   Basis_System.QN1_LL_LL_Stored[LL_site + 1].resize(dim_renorm);
   
   for (int i = 0; i < dim_renorm; i++) {
      Basis_System.QN1_LL_LL_Stored[LL_site + 1][i] = T_Mat.QN1[i];
   }
   
}
