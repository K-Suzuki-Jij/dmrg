#include <ios>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Output_Intersite_Values(std::vector<double> &Val, std::string file_name) {
   
   mkdir("./result", 0777);
   mkdir("./result/CorrelationFunctions", 0777);
   
   std::string Out_Name = "./result/CorrelationFunctions/" + file_name;
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
   for (int site = site_cf_ref; site < system_size; site++) {
      file << std::noshowpos << std::left << std::setw(2) << site_cf_ref        << "  ";
      file << std::noshowpos << std::left << std::setw(2) << site               << "  ";
      file << std::noshowpos << std::left << std::setw(2) << site - site_cf_ref << "  ";
      file << std::showpos   << Val[site - site_cf_ref] << "\n";
   }
   
   file << "\n";
   file.close();
   
}
