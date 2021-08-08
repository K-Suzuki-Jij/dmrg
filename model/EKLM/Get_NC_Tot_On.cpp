//
//  Created by Kohei Suzuki on 2021/01/04.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_NC_Tot_On(CRS &M, double coeef) {
   
   CRS Temp1, Temp2;
   
   Temp2.row_dim = Find_Dim_Onsite();
   Temp2.col_dim = Find_Dim_Onsite();
   M.row_dim     = Find_Dim_Onsite();
   M.col_dim     = Find_Dim_Onsite();
   
   for (int i = 0; i <= Temp2.row_dim; i++) {
      Temp2.Row.push_back(0);
      M.Row.push_back(0);
   }
   
   for (int i = 0; i < num_ele_orbit; i++) {
      Get_NC_On(Temp1, i, coeef);
      Matrix_Matrix_Sum(Temp1, Temp2, M);
      Temp2 = M;
   }

}

