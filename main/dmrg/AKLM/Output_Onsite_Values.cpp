//
//  Output_Onsite_Values.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/07/10.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Output_Onsite_Values(std::vector<double> &Val, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_AKLM &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << LL_site + 1 << "_1_1_" << RR_site + 1 << "]_" << Dmrg_Param.now_sweep << "/ExpectationValues/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
   file << "###";
   file << "BC="      << Model.BC;
   file << ",N="      << Model.system_size;
   file << ",Nc="     << Model.tot_ele;
   file << ",Spin="   << Model.lspin*0.5;
   file << ",Sz="     << Model.tot_sz*0.5;
   file << ",Copy="   << Dmrg_Param.Enviro_Copy;
   file << ",Guess="  << Dmrg_Param.Initial_Guess;
   file << ",sweep="  << Dmrg_Param.tot_sweep;
   file << ",Cutoff=" << Dmrg_Param.max_dim_system;
   file << "\n";
   
   file << "###";
   file << "t="     << Model.t;
   file << ",J_xy=" << Model.J_xy;
   file << ",J_z="  << Model.J_z;
   file << ",I_xy=" << Model.I_xy;
   file << ",I_z="  << Model.I_z;
   file << ",K_xy=" << Model.K_xy;
   file << ",K_z="  << Model.K_z;
   file << ",U="    << Model.U;
   file << ",V="    << Model.V;
   file << ",D_z="  << Model.D_z;
   file << ",hc_z=" << Model.hc_z;
   file << ",hl_z=" << Model.hl_z;
   file << ",mu="   << Model.mu;
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
