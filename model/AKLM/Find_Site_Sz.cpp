#include <iostream>
#include "Model_1D_AKLM.hpp"

int Model_1D_AKLM::Find_Site_Sz(int basis) {
   
   ///////////////////////////////////////
   // # <-> [Cherge  ] -- (N,  2*sz)
   // 0 <-> [        ] -- (0,  0   )
   // 1 <-> [up      ] -- (1,  1   )
   // 2 <-> [down    ] -- (1, -1   )
   // 3 <-> [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int basis_lspin = Find_Basis_Lspin(basis);
   int basis_ele   = Find_Basis_Ele(basis);
   int sz_lspin    = lspin - 2*basis_lspin;
   
   if (basis_ele == 0 || basis_ele == 3) {
      return sz_lspin;
   }
   else if (basis_ele == 1) {
      return sz_lspin + 1;
   }
   else if (basis_ele == 2) {
      return sz_lspin - 1;
   }
   else {
      std::cout << "Error in Find_Site_Sz" << std::endl;
      exit(1);
   }
   
}
