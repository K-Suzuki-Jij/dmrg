//
//  Make_Elem_Ham_LLLR.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_LLLR(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLR, Block_Operator &System, std::vector<int> &Ele_LL, Model_1D_AKLM_TVF &Model) {
   
   int LL_site = A_Basis.LL_site;
   std::vector<int> Dummy;//This wont be used
   
   //Onsite Ham
   DMRG_Make_Elem_LL_LLLR(Basis, A_Basis, Inv_LLLR, System.Ham[LL_site], 1.0);
   DMRG_Make_Elem_LR_LLLR(Basis, A_Basis, Inv_LLLR, Model.Ham_On       , 1.0);
   
   //Interaction LL-LR
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.CEven_D_RE[LL_site], Model.CEven_On  , Model.t       , Ele_LL, "LL_LR");
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.CEven_RE[LL_site]  , Model.CEven_D_On, Model.t       , Ele_LL, "LR_LL");
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.COdd_D_RE[LL_site] , Model.COdd_On   , Model.t       , Ele_LL, "LL_LR");
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.COdd_RE[LL_site]   , Model.COdd_D_On , Model.t       , Ele_LL, "LR_LL");
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.SzL_RE[LL_site]    , Model.SzL_On    , Model.I_z     , Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.SpL_RE[LL_site]    , Model.SmL_On    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.SmL_RE[LL_site]    , Model.SpL_On    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.SzC_RE[LL_site]    , Model.SzC_On    , Model.K_z     , Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.SpC_RE[LL_site]    , Model.SmC_On    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.SmC_RE[LL_site]    , Model.SpC_On    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.NC_RE[LL_site]     , Model.NC_On     , Model.V       , Dummy , "No"   );

   DMRG_Clear_Check_Basis(A_Basis, Inv_LLLR, "LLLR");
   
}
