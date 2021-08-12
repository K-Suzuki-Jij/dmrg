//
//  Make_Elem_Ham_LRRL.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/11/16.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_LRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LRRL, Model_1D_XXZ &Model) {
   
   std::vector<int> Dummy;//This wont be used
   
   //Interaction LR-RL
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.Sm_On, Model.Sp_On, 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.Sp_On, Model.Sm_On, 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.Sz_On, Model.Sz_On,     Model.J_z , Dummy, "No");

   DMRG_Clear_Check_Basis(A_Basis, Inv_LRRL);
   
}
