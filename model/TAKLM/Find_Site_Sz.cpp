#include <iostream>
#include "Model_1D_TAKLM.hpp"

int Model_1D_TAKLM::Find_Site_Sz(int basis) {
   
   ///////////////////////////////////////
   // # <-> [Cherge  ] -- (N,  2*sz)
   // 0 <-> [        ] -- (0,  0   )
   // 1 <-> [up      ] -- (1,  1   )
   // 2 <-> [down    ] -- (1, -1   )
   // 3 <-> [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int basis_lspin = Find_Basis_Lspin(basis);
   int basis_ele_1 = Find_Basis_Ele_1(basis);
   int basis_ele_2 = Find_Basis_Ele_2(basis);
   int sz_lspin    = lspin - 2*basis_lspin;
   int sz_ele_1, sz_ele_2;
   
   if (basis_ele_1 == 0 || basis_ele_1 == 3) {
      sz_ele_1 = 0;
   }
   else if (basis_ele_1 == 1) {
      sz_ele_1 = 1;
   }
   else if (basis_ele_1 == 2) {
      sz_ele_1 = -1;
   }
   else {
      std::cout << "Error in Find_Site_Sz" << std::endl;
      exit(1);
   }
   
   if (basis_ele_2 == 0 || basis_ele_2 == 3) {
      sz_ele_2 = 0;
   }
   else if (basis_ele_2 == 1) {
      sz_ele_2 = 1;
   }
   else if (basis_ele_2 == 2) {
      sz_ele_2 = -1;
   }
   else {
      std::cout << "Error in Find_Site_Sz" << std::endl;
      exit(1);
   }
   
   return sz_ele_1 + sz_ele_2 + sz_lspin;
   
}
