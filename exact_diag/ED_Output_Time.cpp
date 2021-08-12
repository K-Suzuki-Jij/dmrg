#include <ios>
#include <fstream>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include "ED.hpp"

void ED_Output_Time(ED_Time &Time, std::string file_name) {
   
   mkdir("./result", 0777);
   std::string Out_Name = "./result/" + file_name;
   std::ofstream file(Out_Name, std::ios::app);
   
   file << std::fixed << std::setprecision(2);
   file << "Bases: " << Time.bases         << "[sec]" << "\n";
   file << "Ham  : " << Time.ham           << "[sec]" << "\n";
   file << "Diag : " << Time.diag          << "[sec]" << "\n";
   file << "Exp  : " << Time.onsite_exp    << "[sec]" << "\n";
   file << "CF   : " << Time.intersite_exp << "[sec]" << "\n";
   file << "Total: " << Time.total         << "[sec]" << "\n";
   file << "\n";
   file.close();
   
}
