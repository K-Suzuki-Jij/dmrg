//
//  DMRG.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/01.
//

#include "Header.hpp"

void DMRG(Model_1D_EKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param) {

   Model.Get_Onsite_Operator();

   Block_Operator System, Enviro;
   Allocate_Block(System, Model);
   Allocate_Block(Enviro, Model);
   
   //Print_CRS(Model.CUp_On[0], "ddd");
   //std::exit(0);
   
   DMRG_Basis_Stored Basis_System, Basis_Enviro;
   Allocate_Basis(Basis_System, Model);
   Allocate_Basis(Basis_Enviro, Model);
   
   Dmrg_Param.renorm_now_iter = 0;
   Dmrg_Param.renorm_tot_iter = Model.system_size/2 - 1 + Dmrg_Param.tot_sweep*(Model.system_size - 4);
   Dmrg_Param.now_sweep = 0;
   
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
