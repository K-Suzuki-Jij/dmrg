//
//  Renormalize_System.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Renormalize_System(Block_Operator &System, DMRG_LLLR_Basis &Basis_LLLR, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_System, DMRG_T_Mat &T_Mat, Model_1D_AKLM_TVF &Model) {
   
   int LL_site = Block_Compo.LL_site;
   
   CRS Ham_LLLR;
   Get_Ham_LLLR(System, Basis_LLLR, Block_Compo, Basis_System.QN1_LL_LL_Stored[LL_site], Model, Ham_LLLR);
   
   Change_Matrix_Basis(Ham_LLLR, T_Mat.Trans_Mat, System.Ham[LL_site + 1]);
   
   std::vector<int> Dummy;
   DMRG_Renormalize_Matrix_LR(Model.CEven_On  , T_Mat.Trans_Mat, System.CEven_RE[LL_site + 1]  , Basis_System.QN1_LL_LL_Stored[LL_site], "Yes", Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.CEven_D_On, T_Mat.Trans_Mat, System.CEven_D_RE[LL_site + 1], Basis_System.QN1_LL_LL_Stored[LL_site], "Yes", Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.COdd_On   , T_Mat.Trans_Mat, System.COdd_RE[LL_site + 1]   , Basis_System.QN1_LL_LL_Stored[LL_site], "Yes", Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.COdd_D_On , T_Mat.Trans_Mat, System.COdd_D_RE[LL_site + 1] , Basis_System.QN1_LL_LL_Stored[LL_site], "Yes", Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.SzL_On    , T_Mat.Trans_Mat, System.SzL_RE[LL_site + 1]    , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.SpL_On    , T_Mat.Trans_Mat, System.SpL_RE[LL_site + 1]    , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.SmL_On    , T_Mat.Trans_Mat, System.SmL_RE[LL_site + 1]    , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.SzC_On    , T_Mat.Trans_Mat, System.SzC_RE[LL_site + 1]    , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.SpC_On    , T_Mat.Trans_Mat, System.SpC_RE[LL_site + 1]    , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.SmC_On    , T_Mat.Trans_Mat, System.SmC_RE[LL_site + 1]    , Dummy, "No" , Basis_LLLR, Block_Compo);
   DMRG_Renormalize_Matrix_LR(Model.NC_On     , T_Mat.Trans_Mat, System.NC_RE[LL_site + 1]     , Dummy, "No" , Basis_LLLR, Block_Compo);

   int dim_renorm = T_Mat.Trans_Mat.col_dim;
   Basis_System.QN1_LL_LL_Stored[LL_site + 1].resize(dim_renorm);
   Basis_System.QN2_LL_LL_Stored[LL_site + 1].resize(dim_renorm);

   for (int i = 0; i < dim_renorm; i++) {
      Basis_System.QN1_LL_LL_Stored[LL_site + 1][i] = T_Mat.QN1[i];
      Basis_System.QN2_LL_LL_Stored[LL_site + 1][i] = T_Mat.QN2[i];
   }

}
