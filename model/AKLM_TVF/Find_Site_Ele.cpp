#include <iostream>
#include "Model_1D_AKLM_TVF.hpp"

int Model_1D_AKLM_TVF::Find_Site_Ele(int basis) {
   
   ///////////////////////////////////////
   // # <->  [Cherge   ] -- (N,  P)
   // 0 <->  [         ] -- (0,  0)
   // 1 <->  [even     ] -- (1,  0)
   // 2 <->  [odd      ] -- (1,  1)
   // 3 <->  [even, odd] -- (2,  1)
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
