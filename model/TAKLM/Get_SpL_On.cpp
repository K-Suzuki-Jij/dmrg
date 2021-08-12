#include <cmath>
#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SpL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   
   Free_CRS(M);
   M.row_dim = Find_Dim_Onsite();
   M.col_dim = Find_Dim_Onsite();
   
   M.Row.push_back(0);
   
   for (int row_ele_1 = 0; row_ele_1 < dim_ele; row_ele_1++) {
      for (int row_ele_2 = 0; row_ele_2 < dim_ele; row_ele_2++) {
         for (int row_lspin = 0; row_lspin < dim_lspin; row_lspin++) {
            for (int col_ele_1 = 0; col_ele_1 < dim_ele; col_ele_1++) {
               for (int col_ele_2 = 0; col_ele_2 < dim_ele; col_ele_2++) {
                  for (int col_lspin = 0; col_lspin < dim_lspin; col_lspin++) {
                     int a = row_lspin + 1;
                     int b = col_lspin + 1;
                     double val = 0.0;
                     
                     if (a + 1 == b) {
                        val = std::sqrt((lspin*0.5 + 1.0)*(a + b - 1.0) - a*b)*coeef;
                     }
                     
                     if (std::abs(val) > 0.0 && row_ele_1 == col_ele_1 && row_ele_2 == col_ele_2) {
                        M.Val.push_back(val);
                        M.Col.push_back(col_ele_1*dim_ele*dim_lspin + col_ele_2*dim_lspin + col_lspin);
                     }
                  }
               }
            }
            M.Row.push_back(M.Col.size());
         }
      }
   }
      
}
