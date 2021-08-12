//
//  DMRG.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/04.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void DMRG(Model_1D_TAKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param) {
   
   Model.tot_ele_1 = Model.Param_Tot_Ele_1[Dmrg_Param.param_now_iter];
   Model.tot_ele_2 = Model.Param_Tot_Ele_2[Dmrg_Param.param_now_iter];
   Model.tot_sz    = Model.Param_Tot_Sz   [Dmrg_Param.param_now_iter];
   
   Model.t_1    = Model.Param_t_1   [Dmrg_Param.param_now_iter];
   Model.t_2    = Model.Param_t_2   [Dmrg_Param.param_now_iter];
   Model.J_xy_1 = Model.Param_J_xy_1[Dmrg_Param.param_now_iter];
   Model.J_z_1  = Model.Param_J_z_1 [Dmrg_Param.param_now_iter];
   Model.J_xy_2 = Model.Param_J_xy_2[Dmrg_Param.param_now_iter];
   Model.J_z_2  = Model.Param_J_z_2 [Dmrg_Param.param_now_iter];
   Model.I_xy   = Model.Param_I_xy  [Dmrg_Param.param_now_iter];
   Model.I_z    = Model.Param_I_z   [Dmrg_Param.param_now_iter];
   Model.K_xy_1 = Model.Param_K_xy_1[Dmrg_Param.param_now_iter];
   Model.K_z_1  = Model.Param_K_z_1 [Dmrg_Param.param_now_iter];
   Model.K_xy_2 = Model.Param_K_xy_2[Dmrg_Param.param_now_iter];
   Model.K_z_2  = Model.Param_K_z_2 [Dmrg_Param.param_now_iter];
   Model.U_1    = Model.Param_U_1   [Dmrg_Param.param_now_iter];
   Model.V_1    = Model.Param_V_1   [Dmrg_Param.param_now_iter];
   Model.U_2    = Model.Param_U_2   [Dmrg_Param.param_now_iter];
   Model.V_2    = Model.Param_V_2   [Dmrg_Param.param_now_iter];
   Model.D_z    = Model.Param_D_z   [Dmrg_Param.param_now_iter];
   Model.hc_z_1 = Model.Param_hc_z_1[Dmrg_Param.param_now_iter];
   Model.hc_z_2 = Model.Param_hc_z_2[Dmrg_Param.param_now_iter];
   Model.hl_z   = Model.Param_hl_z  [Dmrg_Param.param_now_iter];
   Model.mu_1   = Model.Param_mu_1  [Dmrg_Param.param_now_iter];
   Model.mu_2   = Model.Param_mu_2  [Dmrg_Param.param_now_iter];

   Model.Set_Onsite_Op();
   Model.Check_Parameters();
      
   Block_Operator System, Enviro;
   Allocate_Block(System, Model);
   Allocate_Block(Enviro, Model);
   
   DMRG_Basis_Stored Basis_System, Basis_Enviro;
   Allocate_Basis(Basis_System, Model);
   Allocate_Basis(Basis_Enviro, Model);
   
   Dmrg_Param.renorm_now_iter = 0;
   Dmrg_Param.renorm_tot_iter = Model.system_size/2 - 1 + Dmrg_Param.tot_sweep*(Model.system_size - 4);
   Dmrg_Param.now_sweep = 0;
   
   Print_Status(Dmrg_Param, Model);
   
   DMRG_Block_Information Block;
   DMRG_Basis Basis;
   DMRG_Ground_State GS;
   DMRG_T_Mat T_Mat;
   
   for (Block.LL_site = 0; Block.LL_site < Model.system_size/2 - 1; Block.LL_site++) {
      Block.RR_site = Block.LL_site;
      Renormalize(System, System, Basis_System, Basis_System, Block, Basis, GS, Model, Dmrg_Param, Diag_Param);
   }
   
   for (Block.LL_site = Model.system_size/2 - 1; Block.LL_site < Model.system_size - 4; Block.LL_site++) {
      Block.RR_site = Model.system_size - 4 - Block.LL_site;
      Renormalize(System, System, Basis_System, Basis_System, Block, Basis, GS, Model, Dmrg_Param, Diag_Param);
   }
   
   //Sweep
   if (Dmrg_Param.Enviro_Copy == "Yes") {
      for (Dmrg_Param.now_sweep = 1; Dmrg_Param.now_sweep <= Dmrg_Param.tot_sweep; Dmrg_Param.now_sweep++) {
         for (Block.LL_site = 0; Block.LL_site < Model.system_size - 4; Block.LL_site++) {
            Block.RR_site = Model.system_size - 4 - Block.LL_site;
            Renormalize(System, System, Basis_System, Basis_System, Block, Basis, GS, Model, Dmrg_Param, Diag_Param);
            if (Block.LL_site == Block.RR_site && Dmrg_Param.now_sweep == Dmrg_Param.tot_sweep) {
               break;
            }
         }
      }
   }
   else if (Dmrg_Param.Enviro_Copy == "No") {
      for (Dmrg_Param.now_sweep = 1; Dmrg_Param.now_sweep <= Dmrg_Param.tot_sweep; Dmrg_Param.now_sweep++) {
         for (Block.LL_site = 0; Block.LL_site < Model.system_size - 4; Block.LL_site++) {
            Block.RR_site = Model.system_size - 4 - Block.LL_site;
            if (Block.LL_site == 0) {
               std::swap(System, Enviro);
               std::swap(Basis_System, Basis_Enviro);
            }
            Renormalize(System, Enviro, Basis_System, Basis_Enviro, Block, Basis, GS, Model, Dmrg_Param, Diag_Param);
            if (Block.LL_site == Block.RR_site && Dmrg_Param.now_sweep == Dmrg_Param.tot_sweep) {
               break;
            }
         }
      }
   }
   
}
