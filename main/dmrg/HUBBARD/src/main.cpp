//
//  main.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/10.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

int main(int argc, const char * argv[]) {
   
   Model_1D_HUBBARD Model;
   DMRG_Param    Dmrg_Param;
   Diag_Param    Diag_Param;
   ////////////////////////////////////////////////
   Model.zero_precision = std::pow(10,-15);
   Model.p_thread = 4;
   ////////////////////////////////////////////////
   Model.BC          = "OBC";//"OBC" or "PBC" or "SSD"
   Model.Mat_Type    = "Sym";//"Sym" or "Non_Sym", "Sym" is better than "Non_Sym" but needs large memory
   Model.system_size = 40;
   Model.site_cf_ref = 0;
   ////////////////////////////////////////////////
   Model.Param_Tot_Ele = {20};
   Model.Param_Tot_Sz  = {0};
   ////////////////////////////////////////////////
   Model.Param_t    = {-1.0};
   Model.Param_J_xy = {0.2};
   Model.Param_J_z  = {0.3};
   Model.Param_U    = {0.5};
   Model.Param_V    = {0.2};
   Model.Param_h_z  = {0.0};
   Model.Param_mu   = {0.0};
   ////////////////////////////////////////////////
   Dmrg_Param.Enviro_Copy    = "Yes";
   Dmrg_Param.Initial_Guess  = "Yes";
   Dmrg_Param.tot_sweep      = 3;
   Dmrg_Param.max_dim_system = 800;
   ////////////////////////////////////////////////
   Diag_Param.diag_num      = 1;
   Diag_Param.Diag_Method   = "Lanczos_Slow";
   Diag_Param.diag_acc      = pow(10,-14);
   Diag_Param.diag_min_step = 1;
   Diag_Param.diag_max_step = 800;
   ////////////////////////////////////////////////
   Diag_Param.cg_acc      = pow(10,-6);
   Diag_Param.cg_max_step = 1000;
   ////////////////////////////////////////////////
   Diag_Param.ii_Method   = "CG";
   Diag_Param.ii_acc      = pow(10,-9);
   Diag_Param.ii_diag_add = pow(10,-11);
   Diag_Param.ii_max_step = 3;
   ////////////////////////////////////////////////
   
   Dmrg_Param.param_tot_iter = (int)Model.Param_Tot_Ele.size();
   for (Dmrg_Param.param_now_iter = 0; Dmrg_Param.param_now_iter < Dmrg_Param.param_tot_iter; Dmrg_Param.param_now_iter++) {
      DMRG(Model, Dmrg_Param, Diag_Param);
   }
   
   return 0;
}
