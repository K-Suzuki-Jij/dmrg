//
//  Make_Elem_Ham_LLLRRRRL.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_LLLRRRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, std::vector<int> &Ele_LL, std::vector<int> &Ele_On, std::vector<int> &Ele_RR, Model_1D_AKLM_TVF &Model) {
   
   int LL_site = A_Basis.LL_site;
   int RR_site = A_Basis.RR_site;
   std::vector<int> Dummy;//This wont be used

   //Onsite Ham
   DMRG_Make_Elem_LL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.Ham[LL_site], 1.0);
   DMRG_Make_Elem_LR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Ham_On       , 1.0);
   DMRG_Make_Elem_RR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Enviro.Ham[RR_site], 1.0);
   DMRG_Make_Elem_RL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.Ham_On       , 1.0);

   //Interaction LL-LR
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.CEven_D_RE[LL_site], Model.CEven_On  , Model.t       , Ele_LL, "LL_LR");
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.CEven_RE[LL_site]  , Model.CEven_D_On, Model.t       , Ele_LL, "LR_LL");
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.COdd_D_RE[LL_site] , Model.COdd_On   , Model.t       , Ele_LL, "LL_LR");
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.COdd_RE[LL_site]   , Model.COdd_D_On , Model.t       , Ele_LL, "LR_LL");
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.SzL_RE[LL_site]    , Model.SzL_On    , Model.I_z     , Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.SpL_RE[LL_site]    , Model.SmL_On    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.SmL_RE[LL_site]    , Model.SpL_On    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.SzC_RE[LL_site]    , Model.SzC_On    , Model.K_z     , Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.SpC_RE[LL_site]    , Model.SmC_On    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.SmC_RE[LL_site]    , Model.SpC_On    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_LLLR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, System.NC_RE[LL_site]     , Model.NC_On     , Model.V       , Dummy , "No"   );

   //Interaction RL-RR
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CEven_On  , Enviro.CEven_D_RE[RR_site], Model.t       , Ele_RR, "RR_RL");
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CEven_D_On, Enviro.CEven_RE[RR_site]  , Model.t       , Ele_RR, "RL_RR");
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.COdd_On   , Enviro.COdd_D_RE[RR_site] , Model.t       , Ele_RR, "RR_RL");
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.COdd_D_On , Enviro.COdd_RE[RR_site]   , Model.t       , Ele_RR, "RL_RR");
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SzL_On    , Enviro.SzL_RE[RR_site]    , Model.I_z     , Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SpL_On    , Enviro.SmL_RE[RR_site]    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SmL_On    , Enviro.SpL_RE[RR_site]    , 0.5*Model.I_xy, Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SzC_On    , Enviro.SzC_RE[RR_site]    , Model.K_z     , Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SpC_On    , Enviro.SmC_RE[RR_site]    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SmC_On    , Enviro.SpC_RE[RR_site]    , 0.5*Model.K_xy, Dummy , "No"   );
   DMRG_Make_Elem_RLRR_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.NC_On     , Enviro.NC_RE[RR_site]     , Model.V       , Dummy , "No"   );

   if (A_Basis.c_obc) {
      //Interaction LR-RL
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CEven_D_On, Model.CEven_On  , Model.t       , Ele_RR, Ele_On, "LR_RL");
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.CEven_On  , Model.CEven_D_On, Model.t       , Ele_RR, Ele_On, "RL_LR");
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.COdd_D_On , Model.COdd_On   , Model.t       , Ele_RR, Ele_On, "LR_RL");
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.COdd_On   , Model.COdd_D_On , Model.t       , Ele_RR, Ele_On, "RL_LR");
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SzL_On    , Model.SzL_On    , Model.I_z     , Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SpL_On    , Model.SmL_On    , 0.5*Model.I_xy, Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SmL_On    , Model.SpL_On    , 0.5*Model.I_xy, Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SzC_On    , Model.SzC_On    , Model.K_z     , Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SpC_On    , Model.SmC_On    , 0.5*Model.K_xy, Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.SmC_On    , Model.SpC_On    , 0.5*Model.K_xy, Dummy , Dummy , "No"   );
      DMRG_Make_Elem_LRRL_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL, Model.NC_On     , Model.NC_On     , Model.V       , Dummy , Dummy , "No"   );
   }
   else {
      std::cout << "Error in Make_Elem_Ham_LLLRRRRL" << std::endl;
      exit(1);
   }
   
   DMRG_Make_Elem_Zero_LLLRRRRL(Basis, A_Basis, Inv_LLLRRRRL);
   
   DMRG_Clear_Check_Basis(A_Basis, Inv_LLLRRRRL, "LLLRRRRL");
   
}
