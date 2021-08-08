//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include <iostream>
#include "Model_1D_EKLM.hpp"

int Model_1D_EKLM::Find_Basis_LSpin(int basis, int target_lspin_orbit) {
   
   basis = Find_Basis_LSpin(basis);
   
   for (int lspin_orbit = 0; lspin_orbit < target_lspin_orbit; lspin_orbit++) {
      basis /= (Magnitude_2LSpin[lspin_orbit] + 1);
   }
   
   return basis%(Magnitude_2LSpin[target_lspin_orbit] + 1);
   
}
