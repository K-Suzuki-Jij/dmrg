//
//  Output_Average_Values.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/11.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Output_Average_Values(std::vector<double> &Val, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << LL_site + 1 << "_1_1_" << RR_site + 1 << "]_" << Dmrg_Param.now_sweep << "/AverageValues/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
   file << std::fixed << std::showpos << std::setprecision(5);
   file << Model.t_1    << "  ";
   file << Model.t_2    << "  ";
   file << Model.J_xy_1 << "  ";
   file << Model.J_z_1  << "  ";
   file << Model.J_xy_2 << "  ";
   file << Model.J_z_2  << "  ";
   file << Model.I_xy   << "  ";
   file << Model.I_z    << "  ";
   file << Model.K_xy_1 << "  ";
   file << Model.K_z_1  << "  ";
   file << Model.K_xy_2 << "  ";
   file << Model.K_z_2  << "  ";
   file << Model.U_1    << "  ";
   file << Model.U_2    << "  ";
   file << Model.V_1    << "  ";
   file << Model.V_2    << "  ";
   file << Model.D_z    << "  ";
   file << Model.hc_z_1 << "  ";
   file << Model.hc_z_2 << "  ";
   file << Model.hl_z   << "  ";
   file << Model.mu_1   << "  ";
   file << Model.mu_2   << "  ";
   file << std::setprecision(15);
   file << std::showpos << std::accumulate(Val.begin(), Val.end(), 0.0)/Val.size();
   file << "\n";
   
   
}

void Output_Average_Values(double val, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << LL_site + 1 << "_1_1_" << RR_site + 1 << "]_" << Dmrg_Param.now_sweep << "/AverageValues/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
   file << std::fixed << std::showpos << std::setprecision(5);
   file << Model.t_1    << "  ";
   file << Model.t_2    << "  ";
   file << Model.J_xy_1 << "  ";
   file << Model.J_z_1  << "  ";
   file << Model.J_xy_2 << "  ";
   file << Model.J_z_2  << "  ";
   file << Model.I_xy   << "  ";
   file << Model.I_z    << "  ";
   file << Model.K_xy_1 << "  ";
   file << Model.K_z_1  << "  ";
   file << Model.K_xy_2 << "  ";
   file << Model.K_z_2  << "  ";
   file << Model.U_1    << "  ";
   file << Model.U_2    << "  ";
   file << Model.V_1    << "  ";
   file << Model.V_2    << "  ";
   file << Model.D_z    << "  ";
   file << Model.hc_z_1 << "  ";
   file << Model.hc_z_2 << "  ";
   file << Model.hl_z   << "  ";
   file << Model.mu_1   << "  ";
   file << Model.mu_2   << "  ";
   file << std::setprecision(15);
   file << std::showpos << val/(LL_site + RR_site + 4.0);
   file << "\n";
   
}
