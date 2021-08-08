#include <vector>
#include <cmath>
#include <omp.h>

void Normalize(std::vector<double> &V, int p_threads) {
   
   double norm = 0;
   long   dim  = (long)V.size();
   
#pragma omp parallel for reduction (+: norm) num_threads (p_threads)
   for (long i = 0; i < dim; ++i) {
      norm += V[i]*V[i];
   }
   
   norm = 1.0/std::sqrt(norm);
   
#pragma omp parallel for num_threads (p_threads)
   for (long i = 0; i < dim; ++i) {
      V[i] = norm*V[i];
   }
   
}
