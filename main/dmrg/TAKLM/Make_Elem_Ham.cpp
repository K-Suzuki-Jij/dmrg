//
//  Make_Elem_Ham.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/12/31.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham(std::string Mat_Type, const DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv, const Block_Operator &Block_Op, const std::vector<int> &Ele, const Model_1D_TAKLM &Model) {
   
   int LL_site = A_Basis.LL_site;
   int RR_site = A_Basis.RR_site;
   
   if (Mat_Type == "LLLR") {
      //Onsite Ham
      DMRG_Make_Elem_Onsite_LLLR_RRRL("LL", Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.Ham[LL_site], 1.0);
      DMRG_Make_Elem_Onsite_LLLR_RRRL("LR", Basis.LR, Basis.LL, A_Basis, Inv, Model.Ham_On         , 1.0);

      //Interaction LL-LR
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.CUp_1_D_RE[LL_site]  , Model.CUp_1_On    , Model.t_1       , Ele, "LL_LR");
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.CUp_1_RE[LL_site]    , Model.CUp_1_D_On  , Model.t_1       , Ele, "LR_LL");
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.CDown_1_D_RE[LL_site], Model.CDown_1_On  , Model.t_1       , Ele, "LL_LR");
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.CDown_1_RE[LL_site]  , Model.CDown_1_D_On, Model.t_1       , Ele, "LR_LL");
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.CUp_2_D_RE[LL_site]  , Model.CUp_2_On    , Model.t_2       , Ele, "LL_LR");
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.CUp_2_RE[LL_site]    , Model.CUp_2_D_On  , Model.t_2       , Ele, "LR_LL");
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.CDown_2_D_RE[LL_site], Model.CDown_2_On  , Model.t_2       , Ele, "LL_LR");
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.CDown_2_RE[LL_site]  , Model.CDown_2_D_On, Model.t_2       , Ele, "LR_LL");
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SzL_RE[LL_site]      , Model.SzL_On      , Model.I_z       , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SpL_RE[LL_site]      , Model.SmL_On      , 0.5*Model.I_xy  , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SmL_RE[LL_site]      , Model.SpL_On      , 0.5*Model.I_xy  , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SzC_1_RE[LL_site]    , Model.SzC_1_On    , Model.K_z_1     , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SpC_1_RE[LL_site]    , Model.SmC_1_On    , 0.5*Model.K_xy_1, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SmC_1_RE[LL_site]    , Model.SpC_1_On    , 0.5*Model.K_xy_1, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SzC_2_RE[LL_site]    , Model.SzC_2_On    , Model.K_z_2     , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SpC_2_RE[LL_site]    , Model.SmC_2_On    , 0.5*Model.K_xy_2, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.SmC_2_RE[LL_site]    , Model.SpC_2_On    , 0.5*Model.K_xy_2, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.NC_1_RE[LL_site]     , Model.NC_1_On     , Model.V_1       , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LL, Basis.LR, A_Basis, Inv, Block_Op.NC_2_RE[LL_site]     , Model.NC_2_On     , Model.V_2       , Ele, "No"   );
   }
   else if (Mat_Type == "RRRL") {
      //Onsite Ham
      DMRG_Make_Elem_Onsite_LLLR_RRRL("RR", Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.Ham[RR_site], 1.0);
      DMRG_Make_Elem_Onsite_LLLR_RRRL("RL", Basis.RL, Basis.RR, A_Basis, Inv, Model.Ham_On         , 1.0);

      //Interaction RR-RL
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.CUp_1_D_RE[RR_site]  , Model.CUp_1_On    , Model.t_1       , Ele, "RR_RL");
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.CUp_1_RE[RR_site]    , Model.CUp_1_D_On  , Model.t_1       , Ele, "RL_RR");
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.CDown_1_D_RE[RR_site], Model.CDown_1_On  , Model.t_1       , Ele, "RR_RL");
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.CDown_1_RE[RR_site]  , Model.CDown_1_D_On, Model.t_1       , Ele, "RL_RR");
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.CUp_2_D_RE[RR_site]  , Model.CUp_2_On    , Model.t_2       , Ele, "RR_RL");
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.CUp_2_RE[RR_site]    , Model.CUp_2_D_On  , Model.t_2       , Ele, "RL_RR");
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.CDown_2_D_RE[RR_site], Model.CDown_2_On  , Model.t_2       , Ele, "RR_RL");
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.CDown_2_RE[RR_site]  , Model.CDown_2_D_On, Model.t_2       , Ele, "RL_RR");
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SzL_RE[RR_site]      , Model.SzL_On      , Model.I_z       , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SpL_RE[RR_site]      , Model.SmL_On      , 0.5*Model.I_xy  , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SmL_RE[RR_site]      , Model.SpL_On      , 0.5*Model.I_xy  , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SzC_1_RE[RR_site]    , Model.SzC_1_On    , Model.K_z_1     , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SpC_1_RE[RR_site]    , Model.SmC_1_On    , 0.5*Model.K_xy_1, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SmC_1_RE[RR_site]    , Model.SpC_1_On    , 0.5*Model.K_xy_1, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SzC_2_RE[RR_site]    , Model.SzC_2_On    , Model.K_z_2     , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SpC_2_RE[RR_site]    , Model.SmC_2_On    , 0.5*Model.K_xy_2, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.SmC_2_RE[RR_site]    , Model.SpC_2_On    , 0.5*Model.K_xy_2, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.NC_1_RE[RR_site]     , Model.NC_1_On     , Model.V_1       , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.RR, Basis.RL, A_Basis, Inv, Block_Op.NC_2_RE[RR_site]     , Model.NC_2_On     , Model.V_2       , Ele, "No"   );
   }
   else if (Mat_Type == "LRRL") {
      //Interaction LR-RL
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.CUp_1_D_On  , Model.CUp_1_On    , Model.t_1       , Ele, "LR_RL");
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.CUp_1_On    , Model.CUp_1_D_On  , Model.t_1       , Ele, "RL_LR");
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.CDown_1_D_On, Model.CDown_1_On  , Model.t_1       , Ele, "LR_RL");
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.CDown_1_On  , Model.CDown_1_D_On, Model.t_1       , Ele, "RL_LR");
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.CUp_2_D_On  , Model.CUp_2_On    , Model.t_2       , Ele, "LR_RL");
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.CUp_2_On    , Model.CUp_2_D_On  , Model.t_2       , Ele, "RL_LR");
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.CDown_2_D_On, Model.CDown_2_On  , Model.t_2       , Ele, "LR_RL");
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.CDown_2_On  , Model.CDown_2_D_On, Model.t_2       , Ele, "RL_LR");
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SzL_On      , Model.SzL_On      , Model.I_z       , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SpL_On      , Model.SmL_On      , 0.5*Model.I_xy  , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SmL_On      , Model.SpL_On      , 0.5*Model.I_xy  , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SzC_1_On    , Model.SzC_1_On    , Model.K_z_1     , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SpC_1_On    , Model.SmC_1_On    , 0.5*Model.K_xy_1, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SmC_1_On    , Model.SpC_1_On    , 0.5*Model.K_xy_1, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SzC_2_On    , Model.SzC_2_On    , Model.K_z_2     , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SpC_2_On    , Model.SmC_2_On    , 0.5*Model.K_xy_2, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.SmC_2_On    , Model.SpC_2_On    , 0.5*Model.K_xy_2, Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.NC_1_On     , Model.NC_1_On     , Model.V_1       , Ele, "No"   );
      DMRG_Make_Elem_Intersite_Block(Basis.LR, Basis.RL, A_Basis, Inv, Model.NC_2_On     , Model.NC_2_On     , Model.V_2       , Ele, "No"   );
   }
   else {
      std::cout << "Error in Make_Elem_Ham" << std::endl;
      std::exit(0);
   }
   
   DMRG_Clear_Check_Basis(A_Basis, Inv);
   
}
