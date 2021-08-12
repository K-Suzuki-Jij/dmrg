//
//  Make_Elem_Ham_LLLRRRRL.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_LLLRRRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, std::vector<int> &Ele_LL, std::vector<int> &Ele_On, std::vector<int> &Ele_RR, Model_1D_HUBBARD &Model) {
   
   int LL_site = A_Basis.LL_site;
   int RR_site = A_Basis.RR_site;
   std::vector<int> Dummy;//This wont be used

   //Onsite Ham
   DMRG_Make_Elem_LL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.Ham[LL_site], 1.0);
   DMRG_Make_Elem_LR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Ham_On       , 1.0);
   DMRG_Make_Elem_RR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Enviro.Ham[RR_site], 1.0);
   DMRG_Make_Elem_RL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Ham_On       , 1.0);

   //Interaction LL-LR
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.CUp_D_RE[LL_site]  , Model.CUp_On    , Model.t       , Ele_LL, "LL_LR");
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.CUp_RE[LL_site]    , Model.CUp_D_On  , Model.t       , Ele_LL, "LR_LL");
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.CDown_D_RE[LL_site], Model.CDown_On  , Model.t       , Ele_LL, "LL_LR");
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.CDown_RE[LL_site]  , Model.CDown_D_On, Model.t       , Ele_LL, "LR_LL");
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.Sz_RE[LL_site]     , Model.Sz_On     , Model.J_z     , Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.Sp_RE[LL_site]     , Model.Sm_On     , 0.5*Model.J_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.Sm_RE[LL_site]     , Model.Sp_On     , 0.5*Model.J_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.NC_RE[LL_site]     , Model.NC_On     , Model.V       , Dummy , "No"   );

   //Interaction RL-RR
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CUp_On    , Enviro.CUp_D_RE[RR_site]  , Model.t       , Ele_RR, "RR_RL");
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CUp_D_On  , Enviro.CUp_RE[RR_site]    , Model.t       , Ele_RR, "RL_RR");
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CDown_On  , Enviro.CDown_D_RE[RR_site], Model.t       , Ele_RR, "RR_RL");
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CDown_D_On, Enviro.CDown_RE[RR_site]  , Model.t       , Ele_RR, "RL_RR");
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Sz_On     , Enviro.Sz_RE[RR_site]     , Model.J_z     , Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Sp_On     , Enviro.Sm_RE[RR_site]     , 0.5*Model.J_xy, Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Sm_On     , Enviro.Sp_RE[RR_site]     , 0.5*Model.J_xy, Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.NC_On     , Enviro.NC_RE[RR_site]     , Model.V       , Dummy , "No"   );

   if (A_Basis.c_obc) {
      //Interaction LR-RL
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CUp_D_On  , Model.CUp_On    , Model.t       , Ele_RR, Ele_On, "LR_RL");
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CUp_On    , Model.CUp_D_On  , Model.t       , Ele_RR, Ele_On, "RL_LR");
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CDown_D_On, Model.CDown_On  , Model.t       , Ele_RR, Ele_On, "LR_RL");
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CDown_On  , Model.CDown_D_On, Model.t       , Ele_RR, Ele_On, "RL_LR");
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Sz_On     , Model.Sz_On     , Model.J_z     , Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Sp_On     , Model.Sm_On     , 0.5*Model.J_xy, Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Sm_On     , Model.Sp_On     , 0.5*Model.J_xy, Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.NC_On     , Model.NC_On     , Model.V       , Dummy , Dummy , "No"   );
   }
   else {
      std::cout << "Error in Make_Elem_Ham_LLLRRRRL" << std::endl;
      exit(1);
   }
   
   DMRG_Make_Elem_Zero_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL);
   
   DMRG_Clear_Check_Basis(A_Basis, Inv_LLLRRRRL, "LLLRRRRL");
   
}
