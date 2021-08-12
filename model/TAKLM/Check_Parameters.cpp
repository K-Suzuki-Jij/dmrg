#include <iostream>
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Check_Parameters() {
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
   if (dim_lspin <= 1) {
      std::cout << "Error in Check_Parameters" << std::endl;
      std::cout << "dim_lspin=" << dim_lspin << std::endl;
      exit(1);
   }
   if (dim_ele != 4) {
      std::cout << "Error in Check_Parameters" << std::endl;
      std::cout << "dim_ele=" << dim_ele << std::endl;
      exit(1);
   }
   if (dim_onsite <= 0 || dim_onsite != dim_lspin*dim_ele*dim_ele) {
      std::cout << "Error in Check_Parameters" << std::endl;
      std::cout << "dim_onsite=" << dim_onsite << std::endl;
      exit(1);
   }
   
}
