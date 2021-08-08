//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_CDown_On(CRS &M, int target_ele_orbit, double coeef) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ]
   // 0 <->  [        ]
   // 1 <->  [up      ]
   // 2 <->  [down    ]
   // 3 <->  [up&down ]
   ///////////////////////////////////////
   
   int dim_onsite = Find_Dim_Onsite();
   
   M.row_dim = dim_onsite;
   M.col_dim = dim_onsite;
   
   M.Row.push_back(0);
   
   for (int row = 0; row < dim_onsite; row++) {
      
      int ele_num = 0;
      for (int ele_orbit = 0; ele_orbit < target_ele_orbit; ele_orbit++) {
         int target_row_ele = Find_Basis_Ele(row, ele_orbit);
         if (target_row_ele == 1 || target_row_ele == 2) {
            ele_num += 1;
         }
         else if (target_row_ele == 3) {
            ele_num += 2;
         }
      }
      
      int sign = 1;
      if (ele_num%2 == 1) {
         sign = -1;
      }
      
      for (int col = 0; col < dim_onsite; col++) {
         int c_other = 1;
         
         for (int ele_orbit = 0; ele_orbit < num_ele_orbit; ele_orbit++) {
            int target_row_ele = Find_Basis_Ele(row, ele_orbit);
            int target_col_ele = Find_Basis_Ele(col, ele_orbit);
            if (ele_orbit != target_ele_orbit && target_row_ele != target_col_ele) {
               c_other = 0;
               break;
            }
         }
       
         int target_row_ele = Find_Basis_Ele(row, target_ele_orbit);
         int target_col_ele = Find_Basis_Ele(col, target_ele_orbit);
         
         int c_ele_1 = (target_row_ele == 0 && target_col_ele == 2);
         int c_ele_2 = (target_row_ele == 1 && target_col_ele == 3);
         
         int sign_in = 1;
         if (c_ele_2) {
            sign_in = -1;
         }
         
         double val = coeef*sign*sign_in;
         
         if (c_other && Find_Basis_LSpin(row) == Find_Basis_LSpin(col) && (c_ele_1 || c_ele_2) && std::abs(val) > zero_precision) {
            M.Val.push_back(val);
            M.Col.push_back(col);
         }
      }
      M.Row.push_back(M.Col.size());
   }
   
  
   
}
