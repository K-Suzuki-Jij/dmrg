//
//  Output_Onsite_Values.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/14.
//

#include "Header.hpp"

void Output_Onsite_Values(const std::vector<double> &Val, std::string file_name, const DMRG_Block_Information &Block, const DMRG_Param &Dmrg_Param, const Model_1D_EKLM &Model) {
   
   std::stringstream Out_Name;
   Out_Name << "./result/[" << Block.LL_site + 1 << "_1_1_" << Block.RR_site + 1 << "]_" << Dmrg_Param.now_sweep << "/ExpectationValues/";
   std::filesystem::create_directories(Out_Name.str());
   
   Out_Name << file_name;
   std::ofstream file(Out_Name.str(), std::ios::app);
   
   file << "###";
   file << "BC="      << Model.BC;
   file << ",N="      << Model.system_size;
   for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
      file << ", Nc_" << ele_orbit << "=" << Model.Tot_Ele[ele_orbit];
   }
   file << ", Sz=" << Model.tot_sz;
   
   for (int spin_orbit = 0; spin_orbit < Model.num_lspin_orbit; spin_orbit++) {
      file << ", Spin_" << spin_orbit << "=" << Model.Magnitude_2LSpin[spin_orbit]*0.5;
   }
   file << ",Copy="   << Dmrg_Param.Enviro_Copy;
   file << ",Guess="  << Dmrg_Param.Initial_Guess;
   file << ",sweep="  << Dmrg_Param.tot_sweep;
   file << ",Cutoff=" << Dmrg_Param.max_dim_system;
   file << std::endl;
   
   file << "###";
   for (auto i = 0; i < Model.Ele_Ele_t_site.size(); i++) {
      file << "t_" << i << "=" << Model.Ele_Ele_t_site[i];
   }
   for (auto i = 0; i < Model.Ele_Ele_V_site.size(); i++) {
      file << ", V_" << i << "=" << Model.Ele_Ele_V_site[i];
   }
   for (auto i = 0; i < Model.LSpin_LSpin_Jz_site.size(); i++) {
      file << ", Jz_" << i << "=" << Model.LSpin_LSpin_Jz_site[i];
   }
   for (auto i = 0; i < Model.LSpin_LSpin_Jxy_site.size(); i++) {
      file << ", Jxy_" << i << "=" << Model.LSpin_LSpin_Jxy_site[i];
   }
   
   file << std::fixed << std::setprecision(15);
   for (auto site = 0; site < Val.size(); site++) {
      file << std::noshowpos << std::left << std::setw(2) << Dmrg_Param.param_now_iter << "  ";
      file << std::noshowpos << std::left << std::setw(3) << site << "  ";
      file << std::showpos   << Val[site] << "\n";
   }
   
   file << "\n";
   file.close();
   
}
