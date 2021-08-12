//
//  Renormalize_System.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/05.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Renormalize_System(Block_Operator &System, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, int LL_site, DMRG_T_Mat &T_Mat, Model_1D_TAKLM &Model) {
      
   int dim_LL  = (int)Basis_System.QN1_LL_LL_Stored[LL_site].size();
   std::vector<int> Ele_LL(dim_LL, 0);
   
   for (int i = 0; i < dim_LL; i++) {
      Ele_LL[i] = Basis_System.QN1_LL_LL_Stored[LL_site][i] + Basis_System.QN2_LL_LL_Stored[LL_site][i];
   }
      
   CRS Ham_LLLR;
   Get_Ham_Block("LLLR", System, Basis, Ele_LL, LL_site, int(), Model, Ham_LLLR);
   Change_Matrix_Basis(Ham_LLLR, T_Mat.Trans_Mat, System.Ham[LL_site + 1], Model.p_thread);
   
   std::vector<int> Dummy;
#pragma omp parallel sections num_threads (Model.p_thread)
   {
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.CUp_1_On    , T_Mat.Trans_Mat, System.CUp_1_RE[LL_site + 1]    , Ele_LL, "Yes", Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.CUp_1_D_On  , T_Mat.Trans_Mat, System.CUp_1_D_RE[LL_site + 1]  , Ele_LL, "Yes", Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.CDown_1_On  , T_Mat.Trans_Mat, System.CDown_1_RE[LL_site + 1]  , Ele_LL, "Yes", Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.CDown_1_D_On, T_Mat.Trans_Mat, System.CDown_1_D_RE[LL_site + 1], Ele_LL, "Yes", Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.CUp_2_On    , T_Mat.Trans_Mat, System.CUp_2_RE[LL_site + 1]    , Ele_LL, "Yes", Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.CUp_2_D_On  , T_Mat.Trans_Mat, System.CUp_2_D_RE[LL_site + 1]  , Ele_LL, "Yes", Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.CDown_2_On  , T_Mat.Trans_Mat, System.CDown_2_RE[LL_site + 1]  , Ele_LL, "Yes", Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.CDown_2_D_On, T_Mat.Trans_Mat, System.CDown_2_D_RE[LL_site + 1], Ele_LL, "Yes", Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SzL_On      , T_Mat.Trans_Mat, System.SzL_RE[LL_site + 1]      , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SpL_On      , T_Mat.Trans_Mat, System.SpL_RE[LL_site + 1]      , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SmL_On      , T_Mat.Trans_Mat, System.SmL_RE[LL_site + 1]      , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SzC_1_On    , T_Mat.Trans_Mat, System.SzC_1_RE[LL_site + 1]    , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SpC_1_On    , T_Mat.Trans_Mat, System.SpC_1_RE[LL_site + 1]    , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SmC_1_On    , T_Mat.Trans_Mat, System.SmC_1_RE[LL_site + 1]    , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.NC_1_On     , T_Mat.Trans_Mat, System.NC_1_RE[LL_site + 1]     , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SzC_2_On    , T_Mat.Trans_Mat, System.SzC_2_RE[LL_site + 1]    , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SpC_2_On    , T_Mat.Trans_Mat, System.SpC_2_RE[LL_site + 1]    , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.SmC_2_On    , T_Mat.Trans_Mat, System.SmC_2_RE[LL_site + 1]    , Dummy , "No" , Basis.LLLR);
#pragma omp section
      DMRG_Renormalize_Matrix_LR(Model.NC_2_On     , T_Mat.Trans_Mat, System.NC_2_RE[LL_site + 1]     , Dummy , "No" , Basis.LLLR);
   }
   
   int dim_renorm = T_Mat.Trans_Mat.col_dim;
   Basis_System.QN1_LL_LL_Stored[LL_site + 1].resize(dim_renorm);
   Basis_System.QN2_LL_LL_Stored[LL_site + 1].resize(dim_renorm);
   Basis_System.QN3_LL_LL_Stored[LL_site + 1].resize(dim_renorm);
   
   for (int i = 0; i < dim_renorm; i++) {
      Basis_System.QN1_LL_LL_Stored[LL_site + 1][i] = T_Mat.QN1[i];
      Basis_System.QN2_LL_LL_Stored[LL_site + 1][i] = T_Mat.QN2[i];
      Basis_System.QN3_LL_LL_Stored[LL_site + 1][i] = T_Mat.QN3[i];
   }   
   
}
