//
//  Renormalize_System.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Renormalize_System(Block_Operator &System, DMRG_LLLR_Basis &Basis_LLLR, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_System, DMRG_T_Mat &T_Mat, Model_1D_HUBBARD &Model) {
   
   int LL_site = Block_Compo.LL_site;
   
   CRS Ham_LLLR;
   Get_Ham_LLLR(System, Basis_LLLR, Block_Compo, Basis_System.QN1_LL_LL_Stored[LL_site], Model, Ham_LLLR);
   
   Change_Matrix_Basis(Ham_LLLR, T_Mat.Trans_Mat, System.Ham[LL_site + 1]);
   
   std::vector<int> Dummy;
   DMRG_Renormalize_Matrix_LR(Model.CUp_On    , T_Mat.Trans_Mat, System.CUp_RE[LL_site + 1]    , Basis_System.QN1_LL_LL_Stored[LL_site], "Yes", Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.CUp_D_On  , T_Mat.Trans_Mat, System.CUp_D_RE[LL_site + 1]  , Basis_System.QN1_LL_LL_Stored[LL_site], "Yes", Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.CDown_On  , T_Mat.Trans_Mat, System.CDown_RE[LL_site + 1]  , Basis_System.QN1_LL_LL_Stored[LL_site], "Yes", Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.CDown_D_On, T_Mat.Trans_Mat, System.CDown_D_RE[LL_site + 1], Basis_System.QN1_LL_LL_Stored[LL_site], "Yes", Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.Sz_On     , T_Mat.Trans_Mat, System.Sz_RE[LL_site + 1]     , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.Sp_On     , T_Mat.Trans_Mat, System.Sp_RE[LL_site + 1]     , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.Sm_On     , T_Mat.Trans_Mat, System.Sm_RE[LL_site + 1]     , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.NC_On     , T_Mat.Trans_Mat, System.NC_RE[LL_site + 1]     , Dummy, "No" , Basis_LLLR, Block_Compo);

   int dim_renorm = T_Mat.Trans_Mat.col_dim;
   Basis_System.QN1_LL_LL_Stored[LL_site + 1].resize(dim_renorm);
   Basis_System.QN2_LL_LL_Stored[LL_site + 1].resize(dim_renorm);

   for (int i = 0; i < dim_renorm; i++) {
      Basis_System.QN1_LL_LL_Stored[LL_site + 1][i] = T_Mat.QN1[i];
      Basis_System.QN2_LL_LL_Stored[LL_site + 1][i] = T_Mat.QN2[i];
   }

}
