#include <cmath>
#include <iostream>
#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_SpL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   
   Clear_CRS(M);
   M.row_dim = Find_Dim_Onsite();
   M.col_dim = Find_Dim_Onsite();
   
   M.Row.push_back(0);
   
   for (int row_ele = 0; row_ele < dim_ele; row_ele++) {
      for (int row_lspin = 0; row_lspin < dim_lspin; row_lspin++) {
         for (int col_ele = 0; col_ele < dim_ele; col_ele++) {
            for (int col_lspin = 0; col_lspin < dim_lspin; col_lspin++) {
               
               int mz_row     = lspin - 2*(row_lspin/dim_parity);
               int parity_row = row_lspin%dim_parity;
               int mz_col     = lspin - 2*(col_lspin/dim_parity);
               int parity_col = col_lspin%dim_parity;
               
               double c_p = std::sqrt(std::fabs(lspin*0.5*(lspin*0.5 + 1.0) - mz_col*0.5*(mz_col*0.5 + 1.0)));
               double c_m = std::sqrt(std::fabs(lspin*0.5*(lspin*0.5 + 1.0) - mz_col*0.5*(mz_col*0.5 - 1.0)));
               double val = 0.0;
               
               if (mz_col > 2) {
                  if (mz_col + 2 == mz_row) {
                     val = 0.5*c_p;
                  }
                  else if (mz_col - 2 == mz_row) {
                     if (parity_col == parity_row) {
                        val = 0.5*c_m;
                     }
                     else {
                        val = -0.5*c_m;
                     }
                  }
                  else {
                     val = 0.0;
                  }
               }
               else if (mz_col == 2) {
                  if (mz_col + 2 == mz_row) {
                     val = 0.5*c_p;
                  }
                  else if (mz_col - 2 == mz_row) {
                     if (parity_col == 0) {
                        val = c_m/std::sqrt(2);
                     }
                     else {
                        val = -c_m/std::sqrt(2);
                     }
                  }
                  else {
                     val = 0.0;
                  }
               }
               
               else if (mz_col == 1) {
                  if (mz_col + 2 == mz_row) {
                     val = 0.5*c_p;
                  }
                  else if (2 - mz_col == mz_row) {
                     if (parity_col == 0) {
                        val = 0.5*c_m;
                     }
                     else {
                        val = -0.5*c_m;
                     }
                  }
                  else {
                     val = 0.0;
                  }
               }
               else if (mz_col == 0) {
                  if (mz_col + 2 == mz_row) {
                     val = c_p/std::sqrt(2);
                  }
                  else {
                     val = 0.0;
                  }
               }
               else {
                  std::cout << "Error in Get_SpL_On" << std::endl;
                  exit(1);
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
