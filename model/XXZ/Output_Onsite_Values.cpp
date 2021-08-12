#include <iostream>
#include <ios>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Output_Onsite_Values(std::vector<double> &Val, std::string file_name) {
   
   if ((int)Val.size() != system_size) {
      std::cout << "Error in Model_1D_XXZ::Output_Onsite_Values" << std::endl;
      exit(1);
   }
   
   mkdir("./result", 0777);
   mkdir("./result/ExpectationValues", 0777);
   
   std::string Out_Name = "./result/ExpectationValues/" + file_name;
   std::ofstream file(Out_Name, std::ios::app);
   
   file << "###1D_XXZ_ED";
   file << ", BC="          << BC;
   file << ", system_size=" << system_size;
   file << ", spin="        << spin*0.5;
   file << ", total_sz="    << tot_sz*0.5;
   file << "\n";
   
   file << "###";
   file << "J_xy="  << J_xy;
   file << ", J_z=" << J_z;
   file << ", D_z=" << D_z;
   file << ", h_z=" << h_z;
   file << "\n";
   
   file << std::fixed << std::setprecision(15);
   for (int site = 0; site < system_size; site++) {
      file << std::noshowpos << std::left << std::setw(2) << site << "  ";
      file << std::showpos   << Val[site] << "\n";
   }
   
   file << "\n";
   file.close();
   
}
