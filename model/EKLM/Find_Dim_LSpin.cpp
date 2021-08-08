//
//  Find_Dim_LSpin.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include <iostream>
#include "Model_1D_EKLM.hpp"

int Model_1D_EKLM::Find_Dim_LSpin() {
   
   int dim = 1;
   
   for (int spin_orbit = 0; spin_orbit < num_lspin_orbit; spin_orbit++) {
      dim *= (Magnitude_2LSpin[spin_orbit] + 1);
   }
   
   if (dim <= 0) {
      std::cout << "Error in Find_Dim_LSpin" << std::endl;
      std::cout << "dim=" << dim << std::endl;
      std::exit(0);
   }
   
   return dim;
   
}
