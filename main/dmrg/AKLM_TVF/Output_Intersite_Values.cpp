//
//  Output_Intersite_Values.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Output_Intersite_Values(std::vector<double> &Val_CF, std::vector<double> &Val_On, std::string file_name, DMRG_Block_Composition &Block_Compo, DMRG_Param &Dmrg_Param, Model_1D_AKLM_TVF &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << Block_Compo.LL_site + 1 << "_1_1_" << Block_Compo.RR_site + 1 << "]_" << Dmrg_Param.now_sweep << "/CorrelationFunctions/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
   file << "###";
   file << "BC="      << Model.BC;
   file << ",N="      << Model.system_size;
   file << ",Nc="     << Model.tot_ele;
   file << ",Spin="   << Model.lspin*0.5;
   file << ",P="      << Model.tot_parity;
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
   file << ",hc_x=" << Model.hc_x;
   file << ",hl_x=" << Model.hl_x;
   file << ",mu="   << Model.mu;
   file << "\n";
   
   file << std::fixed << std::setprecision(15);
   for (int site = Model.site_cf_ref; site < Model.system_size; site++) {
      file << std::setfill('0') << std::right << std::noshowpos << std::setw(3) << Dmrg_Param.param_now_iter << " ";
      file << std::setfill('0') << std::right << std::noshowpos << std::setw(3) << Model.site_cf_ref         << " ";
      file << std::setfill('0') << std::right << std::noshowpos << std::setw(3) << site                      << " ";
      file << std::setfill('0') << std::right << std::noshowpos << std::setw(3) << site - Model.site_cf_ref  << " ";
      file << std::showpos   << Val_CF[site - Model.site_cf_ref] << " ";
      file << std::showpos   << Val_CF[site - Model.site_cf_ref] - Val_On[Model.site_cf_ref]*Val_On[site];
      file << "\n";
   }
   
   file << "\n";
   file.close();
   
}
