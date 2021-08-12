#include <iostream>
#include "Model_1D_HUBBARD.hpp"

int Model_1D_HUBBARD::Find_Site_Sz(int basis) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (basis == 0 || basis == 3) {
      return 0;
   }
   else if (basis == 1) {
      return 1;
   }
   else if (basis == 2) {
      return -1;
   }
   else {
      std::cout << "Error in Find_Site_Sz" << std::endl;
      exit(1);
   }
   
}
