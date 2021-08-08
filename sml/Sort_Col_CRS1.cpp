#include "SML.hpp"

void Sort_Col_CRS(CRS &M, int p_threads) {
   
#pragma omp parallel for schedule(guided) num_threads (p_threads)
   for (int i = 0; i < M.row_dim; i++) {
      Quick_Sort_Vector(M.Col, M.Val, M.Row[i], M.Row[i+1]);
   }
  
}

void Sort_Col_CCS(CCS &M, int p_threads) {
   
#pragma omp parallel for schedule(guided) num_threads (p_threads)
   for (int i = 0; i < M.col_dim; i++) {
      Quick_Sort_Vector(M.Row, M.Val, M.Col[i], M.Col[i+1]);
   }
  
}
