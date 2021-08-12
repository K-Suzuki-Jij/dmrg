//
//  Output_Onsite_Values.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/21.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Output_Onsite_Values(std::vector<double> &Val, std::string file_name, DMRG_Block_Composition &Block_Compo, DMRG_Param &Dmrg_Param, Model_1D_HUBBARD &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << Block_Compo.LL_site << "_1_1_" << Block_Compo.RR_site << "]_" << Dmrg_Param.now_sweep << "/ExpectationValues/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
   file << "###";
   file << "BC="      << Model.BC;
   file << ",N="      << Model.system_size;
   file << ",Nc="     << Model.tot_ele;
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
   file << ",U="    << Model.U;
   file << ",V="    << Model.V;
   file << ",h_z="  << Model.h_z;
   file << ",mu="   << Model.mu;
   file << "\n";
   
   file << std::fixed << std::setprecision(15);
   for (size_t site = 0; site < Val.size(); site++) {
      file << std::noshowpos << std::left << std::setw(2) << site << "  ";
      file << std::showpos   << Val[site] << "\n";
   }
   
   file << "\n";
   file.close();
   
}
