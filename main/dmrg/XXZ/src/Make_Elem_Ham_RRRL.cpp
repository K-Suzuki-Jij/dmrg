//
//  Make_Elem_Ham_RRRL.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/11/16.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_RRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_RRRL, Block_Operator &Enviro, Model_1D_XXZ &Model) {
   
   int RR_site = A_Basis.RR_site;
   std::vector<int> Dummy;//This wont be used
   
   //Onsite Ham
   DMRG_Make_Elem_RR_Onsite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.Ham[RR_site], 1.0);
   DMRG_Make_Elem_RL_Onsite_RRRL(Basis, A_Basis, Inv_RRRL, Model.Ham_On       , 1.0);
   
   //Interaction LL-LR
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.Sp_RE[RR_site], Model.Sm_On, 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.Sm_RE[RR_site], Model.Sp_On, 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.Sz_RE[RR_site], Model.Sz_On,     Model.J_z , Dummy, "No");

   DMRG_Clear_Check_Basis(A_Basis, Inv_RRRL);
   
}
