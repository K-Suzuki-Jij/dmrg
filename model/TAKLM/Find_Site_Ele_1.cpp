#include <iostream>
#include "Model_1D_TAKLM.hpp"

int Model_1D_TAKLM::Find_Site_Ele_1(int basis) {
   
   ///////////////////////////////////////
   // # <-> [Cherge  ] -- (N,  2*sz)
   // 0 <-> [        ] -- (0,  0   )
   // 1 <-> [up      ] -- (1,  1   )
   // 2 <-> [down    ] -- (1, -1   )
   // 3 <-> [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int basis_ele_1 = Find_Basis_Ele_1(basis);
     
   if (basis_ele_1 == 0) {
      return 0;
   }
   else if (basis_ele_1 == 1 || basis_ele_1 == 2) {
      return 1;
   }
   else if (basis_ele_1 == 3) {
      return 2;
   }
   else {
      std::cout << "Error in Find_Site_Ele_1" << std::endl;
      exit(1);
   }

   
}
