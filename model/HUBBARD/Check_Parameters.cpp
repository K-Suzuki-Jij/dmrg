#include <iostream>
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Check_Parameters() {
    
    if (tot_ele%2 != tot_sz%2) {
        std::cout << "Error in Check_Parameters" << std::endl;
        std::cout << "tot_ele=" << tot_ele << ", tot_sz=" << tot_sz << std::endl;
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
