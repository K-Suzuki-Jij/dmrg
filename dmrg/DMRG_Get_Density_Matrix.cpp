#include <iostream>
#include "DMRG.hpp"

void DMRG_Get_Density_Matrix(std::vector<std::vector<double>> &M, int block, const std::vector<int> &Dim_Block, const std::vector<double> &GS_Vec, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Basis_LLLR &Basis_LLLR) {
   
   int dim_dm      = Dim_Block[block];
   int dim_LLLR    = (int)Basis_LLLR.LL.size();
   int base_LRRRRL = Basis_LLLRRRRL.base_LRRRRL;
   int base_RRRL   = Basis_LLLRRRRL.base_RRRL;
   int base_RL     = Basis_LLLRRRRL.base_RL;

   std::vector<int> Initial_Index_LLLR(dim_LLLR + 1, 0);
   
   M.resize(dim_dm);
   for (int i = 0; i < dim_dm; i++) {
      M[i].resize(dim_dm);
   }
   
   for (int i = 0; i < dim_LLLR; i++) {
      Initial_Index_LLLR[i + 1] = Initial_Index_LLLR[i] + Basis_LLLR.Count_Enviro[i];
   }
   
   //Check Point
   if (Initial_Index_LLLR[dim_LLLR] != (int)Basis_LLLRRRRL.LL.size()) {
      std::cout << "Error in DMRG_Get_Density_Matrix at 1" << std::endl;
      exit(1);
   }
   
   int base = 0;
   for (int i = 0; i < block; i++) {
      base += Dim_Block[i];
   }
      
   for (int row = 0; row < dim_dm; ++row) {
      int row_base          = row + base;
      int LL_row            = Basis_LLLR.LL[row_base];
      int LR_row            = Basis_LLLR.LR[row_base];
      int loop_count        = Basis_LLLR.Count_Enviro[row_base];
      int initial_index_row = Initial_Index_LLLR[row_base];

      for (int col = row; col < dim_dm; ++col) {
         int col_base          = col + base;
         int initial_index_col = Initial_Index_LLLR[col_base];

         //Check Point
         if (row_base >= dim_LLLR || col_base >= dim_LLLR) {
            std::cout << "Error in DMRG_Get_Density_Matrix at 2" << std::endl;
            exit(1);
         }
         
         //Check Point
         if (Basis_LLLR.Count_Enviro[row_base] != Basis_LLLR.Count_Enviro[col_base]) {
            std::cout << "Error in DMRG_Get_Density_Matrix at 3" << std::endl;
            exit(1);
         }
         
         int LL_col = Basis_LLLR.LL[col_base];
         int LR_col = Basis_LLLR.LR[col_base];
         double val = 0.0;
         for (int i = 0; i < loop_count; ++i) {
            int RR_row = Basis_LLLRRRRL.RR[initial_index_row + i];
            int RR_col = Basis_LLLRRRRL.RR[initial_index_col + i];
            int RL_row = Basis_LLLRRRRL.RL[initial_index_row + i];
            int RL_col = Basis_LLLRRRRL.RL[initial_index_col + i];
            
            //Check Point
            if (RR_row != RR_col || RL_row != RL_col) {
               std::cout << "Error in DMRG_Get_Density_Matrix at 4" << std::endl;
               exit(1);
            }
            
            long LLLRRRRL_row = (long)LL_row*base_LRRRRL + LR_row*base_RRRL + RR_row*base_RL + RL_row;
            long LLLRRRRL_col = (long)LL_col*base_LRRRRL + LR_col*base_RRRL + RR_col*base_RL + RL_col;
            int inv_row = Basis_LLLRRRRL.Inv[LLLRRRRL_row];
            int inv_col = Basis_LLLRRRRL.Inv[LLLRRRRL_col];
            
            //Check Point
            if (inv_row < 0 || inv_col < 0) {
               std::cout << "Error in DMRG_Get_Density_Matrix at 5" << std::endl;
               exit(1);
            }
            
            val += GS_Vec[inv_row]*GS_Vec[inv_col];
            
         }
         M[row][col] = val;
      }
   }
   
   
   for (int row = 0; row < dim_dm; row++) {
      for (int col = row; col < dim_dm; col++) {
         M[col][row] = M[row][col];
      }
   }
   
}
