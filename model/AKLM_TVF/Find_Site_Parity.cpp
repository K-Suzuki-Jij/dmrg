#include <iostream>
#include "Model_1D_AKLM_TVF.hpp"

int Model_1D_AKLM_TVF::Find_Site_Parity(int basis) {
   
   ///////////////////////////////////////
   // # <->  [Cherge   ] -- (N,  P)
   // 0 <->  [         ] -- (0,  0)
   // 1 <->  [even     ] -- (1,  0)
   // 2 <->  [odd      ] -- (1,  1)
   // 3 <->  [even, odd] -- (2,  1)
   ///////////////////////////////////////
   
   int basis_lspin = Find_Basis_Lspin(basis);
   int basis_ele   = Find_Basis_Ele(basis);
   int parity_ele, parity_lspin;
   
   if (basis_ele == 0 || basis_ele == 1) {
      parity_ele = 0;
   }
   else if (basis_ele == 2 || basis_ele == 3) {
      parity_ele = 1;
   }
   else {
      std::cout << "Error in Find_Site_Parity" << std::endl;
      exit(1);
   }
   
   parity_lspin = basis_lspin%dim_parity;
   
   return (parity_ele + parity_lspin)%2;
   
}
