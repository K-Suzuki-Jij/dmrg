//
//  Print_Status.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/07/12.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Print_Status(DMRG_Param &Dmrg_Param, Model_1D_AKLM &Model) {
   
   std::cout << "###";
   std::cout << "BC="      << Model.BC;
   std::cout << ",N="      << Model.system_size;
   std::cout << ",Nc="     << Model.tot_ele;
   std::cout << ",Spin="   << Model.lspin*0.5;
   std::cout << ",Sz="     << Model.tot_sz*0.5;
   std::cout << ",Copy="   << Dmrg_Param.Enviro_Copy;
   std::cout << ",Guess="  << Dmrg_Param.Initial_Guess;
   std::cout << ",sweep="  << Dmrg_Param.tot_sweep;
   std::cout << ",Cutoff=" << Dmrg_Param.max_dim_system;
   std::cout << "\n";
   
   std::cout << "###";
   std::cout << "t="     << Model.t;
   std::cout << ",J_xy=" << Model.J_xy;
   std::cout << ",J_z="  << Model.J_z;
   std::cout << ",I_xy=" << Model.I_xy;
   std::cout << ",I_z="  << Model.I_z;
   std::cout << ",K_xy=" << Model.K_xy;
   std::cout << ",K_z="  << Model.K_z;
   std::cout << ",U="    << Model.U;
   std::cout << ",V="    << Model.V;
   std::cout << ",D_z="  << Model.D_z;
   std::cout << ",hc_z=" << Model.hc_z;
   std::cout << ",hl_z=" << Model.hl_z;
   std::cout << ",mu="   << Model.mu;
   std::cout << "\n";
   
   mkdir("./result", 0777);
   std::string Out_Name = "./result/log.txt";
   std::ofstream file(Out_Name, std::ios::app);
   
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

}
