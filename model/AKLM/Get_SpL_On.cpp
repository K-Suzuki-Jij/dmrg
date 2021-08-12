#include <cmath>
#include "SML.hpp"
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_SpL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   
   Free_CRS(M);
   M.row_dim = Find_Dim_Onsite();
   M.col_dim = Find_Dim_Onsite();
   
   M.Row.push_back(0);
   
   for (int row_ele = 0; row_ele < dim_ele; row_ele++) {
      for (int row_lspin = 0; row_lspin < dim_lspin; row_lspin++) {
         for (int col_ele = 0; col_ele < dim_ele; col_ele++) {
            for (int col_lspin = 0; col_lspin < dim_lspin; col_lspin++) {
               int a = row_lspin + 1;
               int b = col_lspin + 1;
               double val;
               
               if (a + 1 == b) {
                  val = std::sqrt((lspin*0.5 + 1)*(a + b - 1.0) - a*b);
               }
               else {
                  val = 0.0;
               }
               
               val *= coeef;
               
               if (std::fabs(val) > 0.0 && row_ele == col_ele) {
                  M.Val.push_back(val);
                  M.Col.push_back(col_ele*dim_lspin + col_lspin);
               }
            }
         }
         M.Row.push_back(M.Col.size());
      }
   }
}
