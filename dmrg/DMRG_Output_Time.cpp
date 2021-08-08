#include <ios>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "DMRG.hpp"

void DMRG_Output_Time(DMRG_Time &Time) {
   
   mkdir("./result", 0777);
   std::string Out_Name = "./result/time.txt";
   std::ofstream file(Out_Name, std::ios::app);
   
   file << std::fixed;
   file << std::left << std::setw(5) << std::setprecision(1) << Time.total       << "   ";
   file << std::left << std::setw(5) << std::setprecision(1) << Time.diag        << "   ";
   file << std::left << std::setw(5) << std::setprecision(1) << Time.inv_iter    << "   ";
   file << std::left << std::setw(5) << std::setprecision(1) << Time.make_ham    << "   ";
   file << std::left << std::setw(5) << std::setprecision(1) << Time.make_basis  << "   ";
   file << std::left << std::setw(5) << std::setprecision(1) << Time.make_dm_mat << "   ";
   file << std::left << std::setw(5) << std::setprecision(1) << Time.total - (Time.diag + Time.inv_iter + Time.make_ham + Time.make_basis + Time.make_dm_mat);
   file << std::endl;
   
   file.close();
   
   
}
