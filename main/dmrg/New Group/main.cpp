//
//  main.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/02.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

int main(int argc, const char * argv[]) {
   
   Model_1D_TAKLM Model;
   DMRG_Param    Dmrg_Param;
   Diag_Param    Diag_Param;
   ////////////////////////////////////////////////
   Model.zero_precision = std::pow(10,-15);
   Model.p_thread = 8;
   ////////////////////////////////////////////////
   Model.BC          = "OBC";//"OBC" or "PBC" or "SSD"
   Model.system_size = 10;
   Model.lspin       = 1;
   Model.site_cf_ref = 0;
   ////////////////////////////////////////////////
   Model.Param_Tot_Ele_1 = {5};
   Model.Param_Tot_Ele_2 = {5};
   Model.Param_Tot_Sz    = {0};
   ////////////////////////////////////////////////
   Model.Param_t_1    = {-1.0};
   Model.Param_t_2    = {-1.0};
   Model.Param_J_xy_1 = {1.7};
   Model.Param_J_z_1  = {1.7};
   Model.Param_J_xy_2 = {1.7};
   Model.Param_J_z_2  = {1.7};
   Model.Param_I_xy   = {0.5};
   Model.Param_I_z    = {0.6};
   Model.Param_K_xy_1 = {0.7};
   Model.Param_K_z_1  = {0.8};
   Model.Param_K_xy_2 = {0.9};
   Model.Param_K_z_2  = {1.0};
   Model.Param_U_1    = {2.1};
   Model.Param_V_1    = {2.2};
   Model.Param_U_2    = {2.3};
   Model.Param_V_2    = {2.4};
   Model.Param_D_z    = {2.5};
   Model.Param_hc_z_1 = {2.6};
   Model.Param_hc_z_2 = {2.7};
   Model.Param_hl_z   = {2.8};
   Model.Param_mu_1   = {2.9};
   Model.Param_mu_2   = {3.0};
   ////////////////////////////////////////////////
   Dmrg_Param.Enviro_Copy    = "Yes";
   Dmrg_Param.Initial_Guess  = "Yes";
   Dmrg_Param.tot_sweep      = 3;
   Dmrg_Param.max_dim_system = 1200;
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
      
   Dmrg_Param.param_tot_iter = (int)Model.Param_Tot_Ele_1.size();
   for (Dmrg_Param.param_now_iter = 0; Dmrg_Param.param_now_iter < Dmrg_Param.param_tot_iter; Dmrg_Param.param_now_iter++) {
      DMRG(Model, Dmrg_Param, Diag_Param);
   }

   return 0;
}
