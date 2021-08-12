//
//  DMRG.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void DMRG(Model_1D_HUBBARD &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param) {
   
   Model.tot_ele = Model.Param_Tot_Ele[Dmrg_Param.param_now_iter];
   Model.tot_sz  = Model.Param_Tot_Sz [Dmrg_Param.param_now_iter];
   
   Model.t    = Model.Param_t   [Dmrg_Param.param_now_iter];
   Model.J_xy = Model.Param_J_xy[Dmrg_Param.param_now_iter];
   Model.J_z  = Model.Param_J_z [Dmrg_Param.param_now_iter];
   Model.U    = Model.Param_U   [Dmrg_Param.param_now_iter];
   Model.V    = Model.Param_V   [Dmrg_Param.param_now_iter];
   Model.h_z  = Model.Param_h_z[Dmrg_Param.param_now_iter];
   Model.mu   = Model.Param_mu  [Dmrg_Param.param_now_iter];
   
   Model.Check_Parameters();
   Model.Set_Onsite_Op();
   
   Block_Operator System, Enviro;
   Allocate_Block(System, Model);
   Allocate_Block(Enviro, Model);
   
   DMRG_Basis_Stored Basis_System, Basis_Enviro;
   DMRG_Allocate_Bases(Basis_System, Model.system_size);
   DMRG_Allocate_Bases(Basis_Enviro, Model.system_size);
   
   Dmrg_Param.renorm_now_iter = 0;
   Dmrg_Param.renorm_tot_iter = Model.system_size/2 - 1 + Dmrg_Param.tot_sweep*(Model.system_size - 4);
   Dmrg_Param.now_sweep = 0;
   DMRG_Basic_Information Info;
   
   //Print_Status(Dmrg_Param, Model);
   
   for (Info.Block_Compo.LL_site = 0; Info.Block_Compo.LL_site < Model.system_size/2 - 1; Info.Block_Compo.LL_site++) {
      Info.Block_Compo.RR_site = Info.Block_Compo.LL_site;
      Renormalize(System, System, Info, Basis_System, Basis_System, Model, Dmrg_Param, Diag_Param);
   }
   
   for (Info.Block_Compo.LL_site = Model.system_size/2 - 1; Info.Block_Compo.LL_site < Model.system_size - 4; Info.Block_Compo.LL_site++) {
      Info.Block_Compo.RR_site = Model.system_size - 4 - Info.Block_Compo.LL_site;
      Renormalize(System, System, Info, Basis_System, Basis_System, Model, Dmrg_Param, Diag_Param);
   }
   
   //Sweep
   if (Dmrg_Param.Enviro_Copy == "Yes") {
      for (Dmrg_Param.now_sweep = 1; Dmrg_Param.now_sweep <= Dmrg_Param.tot_sweep; Dmrg_Param.now_sweep++) {
         for (Info.Block_Compo.LL_site = 0; Info.Block_Compo.LL_site < Model.system_size - 4; Info.Block_Compo.LL_site++) {
            Info.Block_Compo.RR_site = Model.system_size - 4 - Info.Block_Compo.LL_site;
            Renormalize(System, System, Info, Basis_System, Basis_System, Model, Dmrg_Param, Diag_Param);
            if (Info.Block_Compo.LL_site == Info.Block_Compo.RR_site && Dmrg_Param.now_sweep == Dmrg_Param.tot_sweep) {
               break;
            }
         }
      }
   }
   else if (Dmrg_Param.Enviro_Copy == "No") {
      for (Dmrg_Param.now_sweep = 1; Dmrg_Param.now_sweep <= Dmrg_Param.tot_sweep; Dmrg_Param.now_sweep++) {
         for (Info.Block_Compo.LL_site = 0; Info.Block_Compo.LL_site < Model.system_size - 4; Info.Block_Compo.LL_site++) {
            Info.Block_Compo.RR_site = Model.system_size - 4 - Info.Block_Compo.LL_site;
            if (Info.Block_Compo.LL_site == 0) {
               std::swap(System, Enviro);
               std::swap(Basis_System, Basis_Enviro);
            }
            Renormalize(System, Enviro, Info, Basis_System, Basis_Enviro, Model, Dmrg_Param, Diag_Param);
            if (Info.Block_Compo.LL_site == Info.Block_Compo.RR_site && Dmrg_Param.now_sweep == Dmrg_Param.tot_sweep) {
               break;
            }
         }
      }
   }
   
}
