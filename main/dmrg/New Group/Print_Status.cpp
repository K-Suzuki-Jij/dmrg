//
//  Print_Status.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/04.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Print_Status(DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model) {
   
   std::cout << "###";
   std::cout << "BC="      << Model.BC;
   std::cout << ",N="      << Model.system_size;
   std::cout << ",Nc_1="   << Model.tot_ele_1;
   std::cout << ",Nc_2="   << Model.tot_ele_2;
   std::cout << ",Spin="   << Model.lspin*0.5;
   std::cout << ",Sz="     << Model.tot_sz*0.5;
   std::cout << ",Copy="   << Dmrg_Param.Enviro_Copy;
   std::cout << ",Guess="  << Dmrg_Param.Initial_Guess;
   std::cout << ",sweep="  << Dmrg_Param.tot_sweep;
   std::cout << ",Cutoff=" << Dmrg_Param.max_dim_system;
   std::cout << "\n";
   
   std::cout << "###";
   std::cout << "t_1="     << Model.t_1;
   std::cout << ",t_2="    << Model.t_2;
   std::cout << ",J_xy_1=" << Model.J_xy_1;
   std::cout << ",J_z_1="  << Model.J_z_1;
   std::cout << ",J_xy_2=" << Model.J_xy_2;
   std::cout << ",J_z_2="  << Model.J_z_2;
   std::cout << ",I_xy="   << Model.I_xy;
   std::cout << ",I_z="    << Model.I_z;
   std::cout << ",K_xy_1=" << Model.K_xy_1;
   std::cout << ",K_z_1="  << Model.K_z_1;
   std::cout << ",K_xy_2=" << Model.K_xy_2;
   std::cout << ",K_z_2="  << Model.K_z_2;
   std::cout << ",U_1="    << Model.U_1;
   std::cout << ",U_2="    << Model.U_2;
   std::cout << ",V_1="    << Model.V_1;
   std::cout << ",V_2="    << Model.V_2;
   std::cout << ",D_z="    << Model.D_z;
   std::cout << ",hc_z_1=" << Model.hc_z_1;
   std::cout << ",hc_z_2=" << Model.hc_z_2;
   std::cout << ",hl_z="   << Model.hl_z;
   std::cout << ",mu_1="   << Model.mu_1;
   std::cout << ",mu_2="   << Model.mu_2;
   std::cout << "\n";
   
   mkdir("./result", 0777);
   std::string Out_Name = "./result/log.txt";
   std::ofstream file(Out_Name, std::ios::app);
   
   file << "###";
   file << "BC="      << Model.BC;
   file << ",N="      << Model.system_size;
   file << ",Nc_1="   << Model.tot_ele_1;
   file << ",Nc_2="   << Model.tot_ele_2;
   file << ",Spin="   << Model.lspin*0.5;
   file << ",Sz="     << Model.tot_sz*0.5;
   file << ",Copy="   << Dmrg_Param.Enviro_Copy;
   file << ",Guess="  << Dmrg_Param.Initial_Guess;
   file << ",sweep="  << Dmrg_Param.tot_sweep;
   file << ",Cutoff=" << Dmrg_Param.max_dim_system;
   file << "\n";
   
   file << "###";
   file << "t_1="     << Model.t_1;
   file << ",t_2="    << Model.t_2;
   file << ",J_xy_1=" << Model.J_xy_1;
   file << ",J_z_1="  << Model.J_z_1;
   file << ",J_xy_2=" << Model.J_xy_2;
   file << ",J_z_2="  << Model.J_z_2;
   file << ",I_xy="   << Model.I_xy;
   file << ",I_z="    << Model.I_z;
   file << ",K_xy_1=" << Model.K_xy_1;
   file << ",K_z_1="  << Model.K_z_1;
   file << ",K_xy_2=" << Model.K_xy_2;
   file << ",K_z_2="  << Model.K_z_2;
   file << ",U_1="    << Model.U_1;
   file << ",U_2="    << Model.U_2;
   file << ",V_1="    << Model.V_1;
   file << ",V_2="    << Model.V_2;
   file << ",D_z="    << Model.D_z;
   file << ",hc_z_1=" << Model.hc_z_1;
   file << ",hc_z_2=" << Model.hc_z_2;
   file << ",hl_z="   << Model.hl_z;
   file << ",mu_1="   << Model.mu_1;
   file << ",mu_2="   << Model.mu_2;
   file << "\n";

}
