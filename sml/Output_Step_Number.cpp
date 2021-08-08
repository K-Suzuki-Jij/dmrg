#include <ios>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

void Output_Step_Number(int step_num, double time, std::string file_name) {
   
   mkdir("./SML_log", 0777);
   std::string Out_Name = "./SML_log/" + file_name;
   std::ofstream file(Out_Name, std::ios::app);
   
   file << std::left  << std::setw(4) << step_num << "  ";
   file << std::fixed << std::setprecision(1);
   file << time << "\n";
   
   file.close();
   
}
