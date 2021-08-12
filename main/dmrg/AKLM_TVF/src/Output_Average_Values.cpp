//
//  Output_Average_Values.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/05.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Output_Average_Values(std::vector<double> &Val, std::string file_name, DMRG_Block_Composition &Block_Compo, DMRG_Param &Dmrg_Param, Model_1D_AKLM_TVF &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << Block_Compo.LL_site + 1 << "_1_1_" << Block_Compo.RR_site + 1 << "]_" << Dmrg_Param.now_sweep << "/AverageValues/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
   file << std::fixed << std::showpos << std::setprecision(5);
   file << Model.t    << " ";
   file << Model.J_xy << " ";
   file << Model.J_z  << " ";
   file << Model.I_xy << " ";
   file << Model.I_z  << " ";
   file << Model.K_xy << " ";
   file << Model.K_z  << " ";
   file << Model.U    << " ";
   file << Model.V    << " ";
   file << Model.D_z  << " ";
   file << Model.hc_x << " ";
   file << Model.hl_x << " ";
   file << Model.mu   << " ";
   file << std::setprecision(15);
   file << std::showpos << std::accumulate(Val.begin(), Val.end(), 0.0)/Val.size();
   file << "\n";
   
}

void Output_Average_Values(double val, std::string file_name, DMRG_Block_Composition &Block_Compo, DMRG_Param &Dmrg_Param, Model_1D_AKLM_TVF &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << Block_Compo.LL_site + 1 << "_1_1_" << Block_Compo.RR_site + 1 << "]_" << Dmrg_Param.now_sweep << "/AverageValues/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
   file << std::fixed << std::showpos << std::setprecision(5);
   file << Model.t    << " ";
   file << Model.J_xy << " ";
   file << Model.J_z  << " ";
   file << Model.I_xy << " ";
   file << Model.I_z  << " ";
   file << Model.K_xy << " ";
   file << Model.K_z  << " ";
   file << Model.U    << " ";
   file << Model.V    << " ";
   file << Model.D_z  << " ";
   file << Model.hc_x << " ";
   file << Model.hl_x << " ";
   file << Model.mu   << " ";
   file << std::setprecision(15);
   file << std::showpos << val/(Block_Compo.LL_site + Block_Compo.RR_site + 4.0);
   file << "\n";
   
}
