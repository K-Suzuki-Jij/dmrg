#include <cmath>
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Get_SxSx_On(CRS &M, double coeef) {
   
   Check_Parameters();
   
   Free_CRS(M);
   
   M.row_dim = Find_Dim_Onsite();
   M.col_dim = Find_Dim_Onsite();
   
   M.Row.push_back(0);
   
   for (int row = 0; row < Find_Dim_Onsite(); row++) {
      for (int col = 0; col < Find_Dim_Onsite(); col++) {
         int a = row + 1;
         int b = col + 1;
         double val1,val2,val3,val4;
         
         if (a == b + 2) {
            val1 = std::sqrt((spin*0.5 + 1.0)*(2.0*a - 2.0) - a*(a - 1.0))*std::sqrt((spin*0.5 + 1.0)*(2.0*a - 4.0) - (a - 1.0)*(a - 2.0));
         }
         else {
            val1 = 0;
         }
         
         if (a == b) {
            val2 = (spin*0.5 + 1.0)*(2.0*a - 2.0) - a*(a - 1.0);
            val3 = (spin*0.5 + 1.0)*(2.0*a) - a*(a + 1.0);
         }
         else {
            val2 = 0;
            val3 = 0;
         }
         
         if (a + 2 == b) {
            val4 = std::sqrt((spin*0.5 + 1.0)*(2.0*a) - a*(a + 1.0))*std::sqrt((spin*0.5 + 1.0)*(2.0*a + 2.0) - (a + 2.0)*(a + 1.0));
         }
         else {
            val4 = 0;
         }
         
         double val = 0.25*(val1 + val2 + val3 + val4);
         
         val = val*coeef;
         
         if (std::abs(val) > 0.0) {
            M.Val.push_back(val);
            M.Col.push_back(col);
         }
      }
      M.Row.push_back(M.Col.size());
   }
}
