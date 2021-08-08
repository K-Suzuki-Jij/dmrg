//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include <iostream>
#include "Model_1D_EKLM.hpp"

int Model_1D_EKLM::Find_Basis_Ele(int basis, int target_ele_orbit) {
   
   basis = Find_Basis_Ele(basis);
   
   for (int ele_orbit = 0; ele_orbit < target_ele_orbit; ele_orbit++) {
      basis /= dim_ele_orbit;
   }
   
   return basis%dim_ele_orbit;
   
}
