//
//  Make_Elem_Ham_LLLR.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_LLLR(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLR, Block_Operator &System, std::vector<int> &Ele_LL, Model_1D_HUBBARD &Model) {
   
   int LL_site = A_Basis.LL_site;
   std::vector<int> Dummy;//This wont be used
   
   //Onsite Ham
   DMRG_Make_Elem_LL_LLLR(Basis, A_Basis, Inv_LLLR, System.Ham[LL_site], 1.0);
   DMRG_Make_Elem_LR_LLLR(Basis, A_Basis, Inv_LLLR, Model.Ham_On       , 1.0);
   
   //Interaction LL-LR
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.CUp_D_RE[LL_site]  , Model.CUp_On,     Model.t       , Ele_LL, "LL_LR");
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.CUp_RE[LL_site]    , Model.CUp_D_On,   Model.t       , Ele_LL, "LR_LL");
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.CDown_D_RE[LL_site], Model.CDown_On,   Model.t       , Ele_LL, "LL_LR");
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.CDown_RE[LL_site]  , Model.CDown_D_On, Model.t       , Ele_LL, "LR_LL");
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.Sz_RE[LL_site]     , Model.Sz_On,      Model.J_z     , Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.Sp_RE[LL_site]     , Model.Sm_On,      0.5*Model.J_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.Sm_RE[LL_site]     , Model.Sp_On,      0.5*Model.J_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLR(Basis, A_Basis, Inv_LLLR, System.NC_RE[LL_site]     , Model.NC_On ,     Model.V       , Dummy , "No"   );

   DMRG_Clear_Check_Basis(A_Basis, Inv_LLLR, "LLLR");
   
}
