//
//  main.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include "Header.hpp"

int main(int argc, const char * argv[]) {
   
   Model_1D_EKLM Model;
   Diag_Param    Diag_Param;
   DMRG_Param    Dmrg_Param;
   Model.BC = "OBC";
   Model.p_threads = 4;
   Model.Magnitude_2LSpin = {1};
   Model.Tot_Ele          = {80};
   Model.tot_sz           = 0;
   Model.system_size      = 80;
   Model.site_cf_ref      = 0;
   Model.zero_precision   = std::pow(10, -15);
   ////////////////////////////////////////////////
   Model.Ele_Ele_t_site        = {-1.0};
   Model.Ele_Ele_V_site        = {+0.0};
   Model.LSpin_LSpin_Jz_site   = {+0.0};
   Model.LSpin_LSpin_Jxy_site  = {+0.0};
   Model.ele_mu                = -0.0;
   Model.lspin_Dz              = -0.0;
   Model.ele_lspin_Jz_orbit    = +2.0;
   Model.ele_lspin_Jxy_orbit   = +2.0;
   Model.ele_ele_V_orbit       = -0.0;
   Model.lspin_lspin_Jxy_orbit = +0.0;
   Model.lspin_lspin_Jz_orbit  = +0.0;
   Model.ele_U = +0.0;
   Model.ele_hz = 0;
   Model.lspin_hz = 0;
   ////////////////////////////////////////////////
   Dmrg_Param.Enviro_Copy    = "Yes";
   Dmrg_Param.Initial_Guess  = "Yes";
   Dmrg_Param.tot_sweep      = 2;
   Dmrg_Param.max_dim_system = 1000;
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
   
   Dmrg_Param.param_tot_iter = 1;
   for (Dmrg_Param.param_now_iter = 0; Dmrg_Param.param_now_iter < Dmrg_Param.param_tot_iter; Dmrg_Param.param_now_iter++) {
      DMRG(Model, Dmrg_Param, Diag_Param);
   }
   
   
   return 0;
}
