//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include <iostream>
#include "Model_1D_EKLM.hpp"

int Model_1D_EKLM::Find_Dim_Ele() {
   
   int dim = 1;
   
   for (int ele_orbit = 0; ele_orbit < num_ele_orbit; ele_orbit++) {
      dim *= dim_ele_orbit;
   }
   
   if (dim <= 0) {
      std::cout << "Error in Find_Dim_Ele" << std::endl;
      std::cout << "dim=" << dim << std::endl;
      std::exit(0);
   }
   
   return dim;
   
}
