//
//  Make_Elem_Ham_LRRL.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/12/12.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_LRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LRRL, std::vector<int> &Ele_LR, Model_1D_AKLM &Model) {
   
   std::vector<int> Dummy;//This wont be used
   
   //Interaction LL-LR
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.CUp_D_On  , Model.CUp_On    , Model.t       , Ele_LR, "LR_RL");
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.CUp_On    , Model.CUp_D_On  , Model.t       , Ele_LR, "RL_LR");
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.CDown_D_On, Model.CDown_On  , Model.t       , Ele_LR, "LR_RL");
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.CDown_On  , Model.CDown_D_On, Model.t       , Ele_LR, "RL_LR");
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.SzL_On    , Model.SzL_On    , Model.I_z     , Dummy , "No"   );
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.SpL_On    , Model.SmL_On    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.SmL_On    , Model.SpL_On    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.SzC_On    , Model.SzC_On    , Model.K_z     , Dummy , "No"   );
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.SpC_On    , Model.SmC_On    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.SmC_On    , Model.SpC_On    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_LRRL_Intersite_LRRL(Basis, A_Basis, Inv_LRRL, Model.NC_On     , Model.NC_On     , Model.V       , Dummy , "No"   );

   DMRG_Clear_Check_Basis(A_Basis, Inv_LRRL);
   
}
