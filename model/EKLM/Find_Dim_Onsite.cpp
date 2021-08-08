//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include <iostream>
#include "Model_1D_EKLM.hpp"

int Model_1D_EKLM::Find_Dim_Onsite() {
   
   int dim = Find_Dim_Ele()*Find_Dim_LSpin();

   if (dim <= 1) {
      std::cout << "Error in Find_Dim_Onsite" << std::endl;
      std::cout << "dim=" << dim << std::endl;
      std::exit(0);
   }
   
   return dim;
   
}
