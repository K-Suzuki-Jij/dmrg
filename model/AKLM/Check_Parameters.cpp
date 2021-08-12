#include <iostream>
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Check_Parameters() {
   if (lspin <= 0) {
      std::cout << "Error in Check_Parameters" << std::endl;
      std::cout << "2*lspin=" << lspin << std::endl;
      exit(1);
   }
   if (system_size <= 0) {
      std::cout << "Error in Check_Parameters" << std::endl;
      std::cout << "system_size=" << system_size << std::endl;
      exit(1);
   }
  if (site_cf_ref < 0) {
     std::cout << "Error in Check_Parameters" << std::endl;
     std::cout << "site_cf_ref=" << site_cf_ref << std::endl;
     exit(1);
  }
   
}
