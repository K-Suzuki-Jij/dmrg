#include "SML.hpp"

void Matrix_Constant_Multiplication(CRS &M, double coeef, int p_threads) {
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < M.row_dim; i++) {
      for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
         M.Val[j] *= coeef;
      }
   }
   
}
