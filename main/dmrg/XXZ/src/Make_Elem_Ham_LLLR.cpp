//
//  Make_Elem_Ham_LLLR.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/06/19.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_LLLR(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLR, Block_Operator &System, Model_1D_XXZ &Model) {
   
   int LL_site = A_Basis.LL_site;
   std::vector<int> Dummy;//This wont be used
   
   //Onsite Ham
   DMRG_Make_Elem_LL_Onsite_LLLR(Basis, A_Basis, Inv_LLLR, System.Ham[LL_site], 1.0);
   DMRG_Make_Elem_LR_Onsite_LLLR(Basis, A_Basis, Inv_LLLR, Model.Ham_On, 1.0);
   
   //Interaction LL-LR
   DMRG_Make_Elem_LLLR_Intersite_LLLR(Basis, A_Basis, Inv_LLLR, System.Sp_RE[LL_site], Model.Sm_On, 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_LLLR_Intersite_LLLR(Basis, A_Basis, Inv_LLLR, System.Sm_RE[LL_site], Model.Sp_On, 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_LLLR_Intersite_LLLR(Basis, A_Basis, Inv_LLLR, System.Sz_RE[LL_site], Model.Sz_On,     Model.J_z , Dummy, "No");
      
   DMRG_Clear_Check_Basis(A_Basis, Inv_LLLR);
   
}
