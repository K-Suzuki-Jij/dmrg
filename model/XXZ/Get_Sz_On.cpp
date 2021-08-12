#include <cmath>
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Get_Sz_On(CRS &M, double coeef) {
   
   Check_Parameters();
   
   Free_CRS(M);
   
   M.row_dim = Find_Dim_Onsite();
   M.col_dim = Find_Dim_Onsite();
   
   M.Row.push_back(0);
   
   for (int row = 0; row < Find_Dim_Onsite(); row++) {
      for (int col = 0; col < Find_Dim_Onsite(); col++) {
         int a = row + 1;
         int b = col + 1;
         double val = (spin*0.5 + 1.0 - b)*Delta_Function(a, b);
         val = val*coeef;
         if (std::abs(val) > 0.0) {
            M.Val.push_back(val);
            M.Col.push_back(col);
         }
      }
      M.Row.push_back(M.Col.size());
   }
}
