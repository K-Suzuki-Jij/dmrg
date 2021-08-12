#include <iostream>
#include "Model_1D_AKLM.hpp"

int Model_1D_AKLM::Find_Site_Ele(int basis) {
   
   ///////////////////////////////////////
   // # <-> [Cherge  ] -- (N,  2*sz)
   // 0 <-> [        ] -- (0,  0   )
   // 1 <-> [up      ] -- (1,  1   )
   // 2 <-> [down    ] -- (1, -1   )
   // 3 <-> [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int basis_ele = Find_Basis_Ele(basis);
   
   if (basis_ele == 0) {
      return 0;
   }
   else if (basis_ele == 1 || basis_ele == 2) {
      return 1;
   }
   else if (basis_ele == 3) {
      return 2;
   }
   else {
      std::cout << "Error in Find_Site_Ele" << std::endl;
      exit(1);
   }
   
}
