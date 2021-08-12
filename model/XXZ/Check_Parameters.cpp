#include <iostream>
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Check_Parameters() {
   if (spin <= 0) {
      std::cout << "Error in Check_Parameters" << std::endl;
      std::cout << "2*spin=" << spin << std::endl;
      exit(1);
   }
   if (system_size <= 0) {
      std::cout << "Error in Check_Parameters" << std::endl;
      std::cout << "system_size=" << system_size << std::endl;
      exit(1);
   }
   if (spin%2 == 0 && std::abs(tot_sz)%2 != 0) {
      std::cout << "Error in Check_Parameters" << std::endl;
      std::cout << "2*spin=" << spin;
      std::cout << ", 2*tot_sz=" << tot_sz << std::endl;
      exit(1);
   }
   if (spin%2 == 1) {
      if (system_size%2 == 0 && std::abs(tot_sz)%2 != 0) {
         std::cout << "Error in Check_Parameters" << std::endl;
         std::cout << "2*spin=" << spin;
         std::cout << ", system_size=" << system_size;
         std::cout << ", 2*tot_sz=" << tot_sz << std::endl;
         exit(1);
      }
      if (system_size%2 == 1 && std::abs(tot_sz)%2 != 1) {
         std::cout << "Error in Check_Parameters" << std::endl;
         std::cout << "2*spin=" << spin;
         std::cout << ", system_size=" << system_size;
         std::cout << ", 2*tot_sz=" << tot_sz << std::endl;
         exit(1);
      }
   }
}
