#include <iostream>
#include <climits>
#include <algorithm>
#include <omp.h>
#include "DMRG.hpp"

void DMRG_Get_Trans_Matrix(DMRG_T_Mat &T_Mat, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Basis_LLLR &Basis_LLLR, DMRG_Ground_State &GS, int num_of_qn, int max_dim_system, int p_threads, DMRG_Time &Time) {
   
   double start = omp_get_wtime();
   
   int dim_LLLR  = (int)Basis_LLLR.LL.size();
   int block_num = 1;
   
   if (dim_LLLR == 0) {
      std::cout << "Error in DMRG_Get_Trans_Matrix at 1" << std::endl;
      std::exit(0);
   }
   
   std::vector<int> Dim_Block;
   
   int temp_dim = 0;
   for (int i = 0; i < dim_LLLR - 1; i++) {
      
      for (int j = 0; j < num_of_qn; j++) {
         int c_qn = (Basis_LLLR.QN_LLLR[j][i] > Basis_LLLR.QN_LLLR[j][i+1]);
         for (int k = 0; k < num_of_qn; k++) {
            if (j != k) {
               c_qn = (c_qn && Basis_LLLR.QN_LLLR[k][i] == Basis_LLLR.QN_LLLR[k][i+1]);
            }
         }
         if (c_qn) {
            std::cout << "Error in DMRG_Get_Trans_Matrix at 2" << std::endl;
            std::exit(0);
         }
      }
      
      int flag = 0;
      for (int j = 0; j < num_of_qn; j++) {
         if (Basis_LLLR.QN_LLLR[j][i] != Basis_LLLR.QN_LLLR[j][i+1]) {
            flag = 1;
            break;
         }
      }
      if (flag == 1) {
         Dim_Block.push_back(i - temp_dim + 1);
         temp_dim = i + 1;
         block_num++;
      }
   }
   
   Dim_Block.push_back(dim_LLLR - temp_dim);
   
   std::vector<int> Block_Sorted(dim_LLLR, 0);
   std::vector<int> Basis_Sorted(dim_LLLR, 0);
   std::vector<std::vector<int>> QN_Block_Sorted(num_of_qn, std::vector<int>(dim_LLLR, 0));
   
   int count = 0;
   for (int block = 0; block < block_num; block++) {
      for (int i = 0; i < Dim_Block[block]; i++) {
         Block_Sorted[count] = block;
         Basis_Sorted[count] = i;
         for (int j = 0; j < num_of_qn; j++) {
            QN_Block_Sorted[j][count] = Basis_LLLR.QN_LLLR[j][count];
         }
         count++;
      }
   }
   
   //Check Point
   if (count != dim_LLLR) {
      std::cout << "Error in DMRG_Get_Trans_Matrix at 4" << std::endl;
      std::exit(0);
   }
      
   std::vector<std::vector<std::vector<double>>> Density_Matrix(block_num);
   std::vector<std::vector<std::vector<double>>> Vector_DM(block_num);
   std::vector<std::vector<double>> Value_DM(block_num);
   

#pragma omp parallel for schedule (guided) num_threads (p_threads)
   for (int block = 0; block < block_num; block++) {
      DMRG_Get_Density_Matrix(Density_Matrix[block], block, Dim_Block, GS.Vector, Basis_LLLRRRRL, Basis_LLLR);
      Lapack_Dsyev(Density_Matrix[block], Vector_DM[block], Value_DM[block]);
   }

   DMRG_Sort_Density_Matrix_Bases(Value_DM, Block_Sorted, Basis_Sorted, QN_Block_Sorted, 0, dim_LLLR);

   GS.sum_tr_val = 0;
   T_Mat.Eig_Val_DM.resize(dim_LLLR);
   
   //Density matrix is positive-semidefinite, but we take abs in case of numerical errors
   for (int i = 0; i < dim_LLLR; i++) {
      T_Mat.Eig_Val_DM[i] = std::abs(Value_DM[Block_Sorted[i]][Basis_Sorted[i]]);
      GS.sum_tr_val      += std::abs(Value_DM[Block_Sorted[i]][Basis_Sorted[i]]);
   }
   
   int dim_renorm = DMRG_Find_Renormalized_Dim(Value_DM, Block_Sorted, Basis_Sorted, max_dim_system);
   
   //Store System Information: Truncation Error
   GS.tr_error = 0;
   for (int i = 0; i < dim_renorm; i++) {
      GS.tr_error += std::abs(Value_DM[Block_Sorted[i]][Basis_Sorted[i]]);
   }
   GS.tr_error = std::abs(GS.sum_tr_val - GS.tr_error);
      
   std::vector<int> Row_Base(block_num + 1, 0);
   
   for (int i = 0; i < block_num; i++) {
      Row_Base[i + 1] = Row_Base[i] + Dim_Block[i];
   }

   long elem_num = 0;
   for (int i = 0; i < dim_renorm; i++) {
      elem_num += Dim_Block[Block_Sorted[i]];
   }
   
   T_Mat.Trans_Mat.Col.resize(dim_renorm + 1);
   T_Mat.Trans_Mat.Row.resize(elem_num);
   T_Mat.Trans_Mat.Val.resize(elem_num);
   T_Mat.Trans_Mat.row_dim = dim_LLLR;
   T_Mat.Trans_Mat.col_dim = dim_renorm;
   T_Mat.QN.resize(num_of_qn);
   
   for (int i = 0; i < num_of_qn; i++) {
      T_Mat.QN[i].resize(dim_renorm);
   }

   for (int i = 0; i < dim_renorm; i++) {
      for (int j = 0; j < num_of_qn; j++) {
         T_Mat.QN[j][i] = QN_Block_Sorted[j][i];
      }
   }
   
   elem_num = 0;
   T_Mat.Trans_Mat.Col[0] = 0;
   for (int i = 0; i < dim_renorm; i++) {
      int block_sorted = Block_Sorted[i];
      int basis_sorted = Basis_Sorted[i];
      for (int j = 0; j < Dim_Block[block_sorted]; j++) {
         T_Mat.Trans_Mat.Val[elem_num] = Vector_DM[block_sorted][basis_sorted][j];
         T_Mat.Trans_Mat.Row[elem_num] = Row_Base[block_sorted] + j;
         elem_num++;
      }
      T_Mat.Trans_Mat.Col[i + 1] = elem_num;
   }
   
   if (elem_num != T_Mat.Trans_Mat.Col[dim_renorm]) {
      std::cout << "Error in DMRG_Get_Trans_Matrix at 5" << std::endl;
      exit(1);
   }
   
   Time.make_dm_mat = omp_get_wtime() - start;
 
}
