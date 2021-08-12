//
//  Output_Onsite_Values.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/11.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Output_Onsite_Values(std::vector<double> &Val, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << LL_site + 1 << "_1_1_" << RR_site + 1 << "]_" << Dmrg_Param.now_sweep << "/ExpectationValues/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
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
   
   file << std::fixed << std::setprecision(15);
   for (size_t site = 0; site < Val.size(); site++) {
      file << std::noshowpos << std::left << std::setw(2) << Dmrg_Param.param_now_iter << "  ";
      file << std::noshowpos << std::left << std::setw(3) << site << "  ";
      file << std::showpos   << Val[site] << "\n";
   }
   
   file << "\n";
   file.close();
   
}
