//
//  Created by Kohei Suzuki on 2021/01/01.
//

#include "Model_1D_EKLM.hpp"

int Model_1D_EKLM::Find_Site_Sz(int basis) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ]
   // 0 <->  [        ]
   // 1 <->  [up      ]
   // 2 <->  [down    ]
   // 3 <->  [up&down ]
   ///////////////////////////////////////
   
   int sz_ele = 0;
   for (int ele_orbit = 0; ele_orbit < num_ele_orbit; ele_orbit++) {
      int row_basis = Find_Basis_Ele(basis, ele_orbit);
      if (row_basis == 1) {
         sz_ele += 1;
      }
      else if (row_basis == 2) {
         sz_ele += -1;
      }
   }
   
   int sz_lspin = 0;
   
   for (int lspin_orbit = 0; lspin_orbit < num_lspin_orbit; lspin_orbit++) {
      sz_lspin += Magnitude_2LSpin[lspin_orbit] - 2*Find_Basis_LSpin(basis, lspin_orbit);
   }
   
   return sz_ele + sz_lspin;

}
