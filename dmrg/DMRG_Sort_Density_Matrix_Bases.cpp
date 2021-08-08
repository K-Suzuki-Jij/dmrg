#include <iostream>
#include <omp.h>
#include "DMRG.hpp"

//Quick sort
void DMRG_Sort_Density_Matrix_Bases(std::vector<std::vector<double>> &Value_DM, std::vector<int> &Block_Sorted, std::vector<int> &Basis_Sorted, std::vector<std::vector<int>> &QN_Block_Sorted, long left, long right) {
   
   if (right - left <= 1) {
      return;
   }

   int    num_of_qn   = (int)QN_Block_Sorted.size();
   long   pivot_index = (left + right)/2;
   double pivot       = Value_DM[Block_Sorted[pivot_index]][Basis_Sorted[pivot_index]];

   std::swap(Block_Sorted[pivot_index], Block_Sorted[right - 1]);
   std::swap(Basis_Sorted[pivot_index], Basis_Sorted[right - 1]);
   
   for (int i = 0; i < num_of_qn; i++) {
      std::swap(QN_Block_Sorted[i][pivot_index], QN_Block_Sorted[i][right - 1]);
   }
   
   long index = left;
   
   for (long i = left; i < right - 1; i++) {
      if (Value_DM[Block_Sorted[i]][Basis_Sorted[i]] > pivot) {
         std::swap(Block_Sorted[index], Block_Sorted[i]);
         std::swap(Basis_Sorted[index], Basis_Sorted[i]);
         for (int j = 0; j < num_of_qn; j++) {
            std::swap(QN_Block_Sorted[j][index], QN_Block_Sorted[j][i]);
         }
         index++;
      }
   }
   
   std::swap(Block_Sorted[index], Block_Sorted[right - 1]);
   std::swap(Basis_Sorted[index], Basis_Sorted[right - 1]);
   
   for (int i = 0; i < num_of_qn; i++) {
      std::swap(QN_Block_Sorted[i][index], QN_Block_Sorted[i][right - 1]);
   }
   
   DMRG_Sort_Density_Matrix_Bases(Value_DM, Block_Sorted, Basis_Sorted, QN_Block_Sorted, left     , index);
   DMRG_Sort_Density_Matrix_Bases(Value_DM, Block_Sorted, Basis_Sorted, QN_Block_Sorted, index + 1, right);

}
