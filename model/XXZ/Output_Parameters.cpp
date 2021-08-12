#include <ios>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Output_Parameters(std::string file_name) {
   
   mkdir("./result", 0777);
   std::string Out_Name = "./result/" + file_name;
   std::ofstream file(Out_Name, std::ios::app);
 
   file << "1D_XXZ_ED";
   file << ", BC="          << BC;
   file << ", system_size=" << system_size;
   file << ", spin="        << spin*0.5;
   file << ", total_sz="    << tot_sz*0.5;
   file << ", dim_target="  << dim_target;
   file << "\n";
   
   file << "J_xy="  << J_xy;
   file << ", J_z=" << J_z;
   file << ", D_z=" << D_z;
   file << ", h_z=" << h_z;
   file << "\n";
   
   file << "site_cf_ref="    << site_cf_ref;
   file << ", dim_target="     << dim_target;
   file << ", zero_precision=" << zero_precision;
   file << "\n";
   file.close();
   
}
