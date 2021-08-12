//
//  main.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/05/18.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

int main(int argc, const char * argv[]) {
   
   Model_1D_XXZ Model;
   DMRG_Param   Dmrg_Param;
   Diag_Param   Diag_Param;
   ////////////////////////////////////////////////
   Model.zero_precision = std::pow(10,-15);
   Model.p_thread = 4;
   ////////////////////////////////////////////////
   Model.BC          = "OBC";//"OBC" or "PBC" or "SSD"
   Model.system_size = 20;
   Model.spin        = 1;
   Model.site_cf_ref = 0;
   ////////////////////////////////////////////////
   Model.Param_Tot_Sz = {0};
   ////////////////////////////////////////////////
   Model.Param_J_xy = {1.0};
   Model.Param_J_z  = {1.0};
   Model.Param_D_z  = {0.0};
   Model.Param_h_z  = {0.0};
   ////////////////////////////////////////////////
   Dmrg_Param.Enviro_Copy    = "Yes";
   Dmrg_Param.Initial_Guess  = "Yes";
   Dmrg_Param.tot_sweep      = 3;
   Dmrg_Param.max_dim_system = 800;
   ////////////////////////////////////////////////
   Diag_Param.diag_num      = 1;
   Diag_Param.Diag_Method   = "Lanczos_Slow";
   Diag_Param.diag_acc      = std::pow(10,-14);
   Diag_Param.diag_min_step = 1;
   Diag_Param.diag_max_step = 800;
   ////////////////////////////////////////////////
   Diag_Param.cg_acc      = std::pow(10,-6);
   Diag_Param.cg_max_step = 1000;
   ////////////////////////////////////////////////
   Diag_Param.ii_Method   = "CG";
   Diag_Param.ii_acc      = std::pow(10,-9);
   Diag_Param.ii_diag_add = std::pow(10,-11);
   Diag_Param.ii_max_step = 3;
   ////////////////////////////////////////////////
   
   Dmrg_Param.param_tot_iter = (int)Model.Param_Tot_Sz.size();
   for (Dmrg_Param.param_now_iter = 0; Dmrg_Param.param_now_iter < Dmrg_Param.param_tot_iter; Dmrg_Param.param_now_iter++) {
      DMRG(Model, Dmrg_Param, Diag_Param);
   }
   
   return 0;
}
