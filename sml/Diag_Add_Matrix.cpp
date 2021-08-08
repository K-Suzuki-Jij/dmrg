#include <iostream>
#include "SML.hpp"

void Diag_Add_Matrix(CRS &M, double add, int p_threads) {
   
   if (M.row_dim != M.col_dim) {
      std::cout << "Error in Diag_Add_Matrix at 1" << std::endl;
      std::exit(0);
   }
   
   long count = 0;
   
#pragma omp parallel for reduction (+:count) num_threads (p_threads)
   for (int i = 0; i < M.row_dim; i++) {
      for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
         if (i == M.Col[j]) {
            M.Val[j] += add;
            count++;
            break;
         }
      }
   }
   
   if (count != M.row_dim) {
      std::cout << "Error in Diag_Add_Matrix at 2" << std::endl;
      std::exit(0);
   }
   
}
