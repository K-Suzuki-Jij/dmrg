#include <iostream>
#include "Model_1D_TAKLM.hpp"

int Model_1D_TAKLM::Find_Site_Ele_2(int basis) {
   
   ///////////////////////////////////////
   // # <-> [Cherge  ] -- (N,  2*sz)
   // 0 <-> [        ] -- (0,  0   )
   // 1 <-> [up      ] -- (1,  1   )
   // 2 <-> [down    ] -- (1, -1   )
   // 3 <-> [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int basis_ele_2 = Find_Basis_Ele_2(basis);
     
   if (basis_ele_2 == 0) {
      return 0;
   }
   else if (basis_ele_2 == 1 || basis_ele_2 == 2) {
      return 1;
   }
   else if (basis_ele_2 == 3) {
      return 2;
   }
   else {
      std::cout << "Error in Find_Site_Ele_2" << std::endl;
      exit(1);
   }

   
}
