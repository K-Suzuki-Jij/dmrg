//
//  Make_Elem_Ham_RRRL.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/12/12.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_RRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_RRRL, Block_Operator &Enviro, std::vector<int> &Ele_RR, Model_1D_AKLM &Model) {
   
   int RR_site = A_Basis.RR_site;
   std::vector<int> Dummy;//This wont be used
   
   //Onsite Ham
   DMRG_Make_Elem_RR_Onsite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.Ham[RR_site], 1.0);
   DMRG_Make_Elem_RL_Onsite_RRRL(Basis, A_Basis, Inv_RRRL, Model.Ham_On       , 1.0);
   
   //Interaction RR-RL
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.CUp_D_RE[RR_site]  , Model.CUp_On    , Model.t       , Ele_RR, "RR_RL");
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.CUp_RE[RR_site]    , Model.CUp_D_On  , Model.t       , Ele_RR, "RL_RR");
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.CDown_D_RE[RR_site], Model.CDown_On  , Model.t       , Ele_RR, "RR_RL");
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.CDown_RE[RR_site]  , Model.CDown_D_On, Model.t       , Ele_RR, "RL_RR");
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.SzL_RE[RR_site]    , Model.SzL_On    , Model.I_z     , Dummy , "No"   );
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.SpL_RE[RR_site]    , Model.SmL_On    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.SmL_RE[RR_site]    , Model.SpL_On    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.SzC_RE[RR_site]    , Model.SzC_On    , Model.K_z     , Dummy , "No"   );
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.SpC_RE[RR_site]    , Model.SmC_On    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.SmC_RE[RR_site]    , Model.SpC_On    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_RRRL_Intersite_RRRL(Basis, A_Basis, Inv_RRRL, Enviro.NC_RE[RR_site]     , Model.NC_On     , Model.V       , Dummy , "No"   );

   DMRG_Clear_Check_Basis(A_Basis, Inv_RRRL);
   
}
