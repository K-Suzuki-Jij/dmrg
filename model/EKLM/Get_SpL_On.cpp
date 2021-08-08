//
//  Created by Kohei Suzuki on 2020/12/23.
//

#include <cmath>
#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_SpL_On(CRS &M, int target_lspin_orbit, double coeef) {

   int dim_onsite = Find_Dim_Onsite();
   
   M.row_dim = dim_onsite;
   M.col_dim = dim_onsite;
   
   M.Row.push_back(0);
   
   for (int row = 0; row < dim_onsite; row++) {
      for (int col = 0; col < dim_onsite; col++) {
         
         int c_other = 1;
         
         for (int lspin_orbit = 0; lspin_orbit < num_lspin_orbit; lspin_orbit++) {
            int target_row_lspin = Find_Basis_LSpin(row, lspin_orbit);
            int target_col_lspin = Find_Basis_LSpin(col, lspin_orbit);
            if (lspin_orbit != target_lspin_orbit && target_row_lspin != target_col_lspin) {
               c_other = 0;
               break;
            }
         }
         
         int    a = Find_Basis_LSpin(row, target_lspin_orbit) + 1;
         int    b = Find_Basis_LSpin(col, target_lspin_orbit) + 1;
         double s = Magnitude_2LSpin[target_lspin_orbit]*0.5;

         double val = std::sqrt((s + 1.0)*(a + b - 1) - a*b)*coeef;
         
         if (c_other && Find_Basis_Ele(row) == Find_Basis_Ele(col) && a + 1 == b && std::abs(val) > zero_precision) {
            M.Val.push_back(val);
            M.Col.push_back(col);
         }
      }
      M.Row.push_back(M.Col.size());
   }
   

}

