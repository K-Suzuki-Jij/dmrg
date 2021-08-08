//
//  Created by Kohei Suzuki on 2021/01/01.
//

#include <iostream>
#include "Model_1D_EKLM.hpp"

int Model_1D_EKLM::Find_Site_Ele(int basis, int target_ele_orbit) {
   
   ///////////////////////////////////////
   // # <-> [Cherge  ] -- (N,  2*sz)
   // 0 <-> [        ] -- (0,  0   )
   // 1 <-> [up      ] -- (1,  1   )
   // 2 <-> [down    ] -- (1, -1   )
   // 3 <-> [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int row_basis = Find_Basis_Ele(basis, target_ele_orbit);
   
   if (row_basis == 0) {
      return 0;
   }
   else if (row_basis == 1 || row_basis == 2) {
      return 1;
   }
   else if (row_basis == 3) {
      return 2;
   }
   else {
      std::cout << "Error in Model_1D_EKLM::Find_Site_Ele" << std::endl;
      std::exit(0);
   }

}
