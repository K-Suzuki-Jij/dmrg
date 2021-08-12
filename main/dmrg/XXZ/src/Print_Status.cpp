//
//  Print_Status.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/11/14.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Print_Status(DMRG_Param &Dmrg_Param, Model_1D_XXZ &Model) {
   
   std::cout << "###";
   std::cout << "BC="      << Model.BC;
   std::cout << ",N="      << Model.system_size;
   std::cout << ",Spin="   << Model.spin*0.5;
   std::cout << ",Sz="     << Model.tot_sz*0.5;
   std::cout << ",Copy="   << Dmrg_Param.Enviro_Copy;
   std::cout << ",Guess="  << Dmrg_Param.Initial_Guess;
   std::cout << ",sweep="  << Dmrg_Param.tot_sweep;
   std::cout << ",Cutoff=" << Dmrg_Param.max_dim_system;
   std::cout << "\n";
   
   std::cout << "###";
   std::cout << "J_xy="    << Model.J_xy;
   std::cout << ",J_z="    << Model.J_z;
   std::cout << ",D_z="    << Model.D_z;
   std::cout << ",h_z="    << Model.h_z;
   std::cout << "\n";
   
   mkdir("./result", 0777);
   std::string Out_Name = "./result/log.txt";
   std::ofstream file(Out_Name, std::ios::app);
   
   file << "###";
   file << "BC="      << Model.BC;
   file << ",N="      << Model.system_size;
   file << ",Spin="   << Model.spin*0.5;
   file << ",Sz="     << Model.tot_sz*0.5;
   file << ",Copy="   << Dmrg_Param.Enviro_Copy;
   file << ",Guess="  << Dmrg_Param.Initial_Guess;
   file << ",sweep="  << Dmrg_Param.tot_sweep;
   file << ",Cutoff=" << Dmrg_Param.max_dim_system;
   file << "\n";
   
   file << "###";
   file << "t_1="     << Model.J_xy;
   file << ",t_2="    << Model.J_z;
   file << ",J_xy_1=" << Model.D_z;
   file << ",J_z_1="  << Model.h_z;
   file << "\n";

}
