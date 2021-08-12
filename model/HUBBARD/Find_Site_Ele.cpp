#include <iostream>
#include "Model_1D_HUBBARD.hpp"

int Model_1D_HUBBARD::Find_Site_Ele(int basis) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (basis == 0) {
      return 0;
   }
   else if (basis == 1 || basis == 2) {
      return 1;
   }
   else if (basis == 3) {
      return 2;
   }
   else {
      std::cout << "Error in Find_Site_Ele" << std::endl;
      exit(1);
   }
   
}
