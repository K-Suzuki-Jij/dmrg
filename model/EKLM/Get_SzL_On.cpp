//
//  Created by Kohei Suzuki on 2020/12/23.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_SzL_On(CRS &M, int target_lspin_orbit, double coeef) {

   int dim_onsite = Find_Dim_Onsite();
   
   M.row_dim = dim_onsite;
   M.col_dim = dim_onsite;
   
   M.Row.push_back(0);
   for (int row = 0; row < dim_onsite; row++) {
      for (int col = 0; col < dim_onsite; col++) {
         int target_basis_row_lspin = Find_Basis_LSpin(row, target_lspin_orbit);
         
         double val = (Magnitude_2LSpin[target_lspin_orbit]*0.5 - target_basis_row_lspin)*coeef;
         
         if (Find_Basis_Ele(row) == Find_Basis_Ele(col) && Find_Basis_LSpin(row) == Find_Basis_LSpin(col) && std::abs(val) > zero_precision) {
            M.Val.push_back(val);
            M.Col.push_back(col);
         }
      }
      M.Row.push_back(M.Col.size());
   }
   
}
