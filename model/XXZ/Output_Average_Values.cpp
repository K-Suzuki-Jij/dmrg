#include <ios>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Output_Average_Values(double val, std::string file_name) {
   
   mkdir("./result", 0777);
   mkdir("./result/AverageValues", 0777);
   
   std::string Out_Name = "./result/AverageValues/" + file_name;
   std::ofstream file(Out_Name, std::ios::app);
      
   file << std::noshowpos << std::left  << std::setw(2) << system_size << "  ";
   file << std::fixed     << std::setprecision(1);
   file << std::noshowpos << spin*0.5   << "  ";
   file << std::showpos   << tot_sz*0.5 << "  ";
   file << std::fixed     << std::setprecision(5);
   file << std::showpos   << J_xy << "  ";
   file << std::showpos   << J_z  << "  ";
   file << std::showpos   << D_z  << "  ";
   file << std::showpos   << h_z  << "  ";
   file << std::fixed     << std::setprecision(15);
   file << std::showpos   << val  << "\n";
   
   file.close();
   
}

void Model_1D_XXZ::Output_Average_Values(std::vector<double> &Val, std::string file_name) {
   
   double avg_val = 0.0;
   
   for (size_t i = 0; i < Val.size(); i++) {
      avg_val += Val[i];
   }
   
   avg_val = avg_val/Val.size();
   
   mkdir("./result", 0777);
   mkdir("./result/AverageValues", 0777);
   
   std::string Out_Name = "./result/AverageValues/" + file_name;
   std::ofstream file(Out_Name, std::ios::app);
      
   file << std::noshowpos << std::left  << std::setw(2) << system_size << "  ";
   file << std::fixed     << std::setprecision(1);
   file << std::noshowpos << spin*0.5   << "  ";
   file << std::showpos   << tot_sz*0.5 << "  ";
   file << std::fixed     << std::setprecision(5);
   file << std::showpos   << J_xy << "  ";
   file << std::showpos   << J_z  << "  ";
   file << std::showpos   << D_z  << "  ";
   file << std::showpos   << h_z  << "  ";
   file << std::fixed     << std::setprecision(15);
   file << std::showpos   << avg_val << "\n";
   
   file.close();
   
}
