#include <cmath>
#include <vector>
#include <omp.h>

double L2_Norm(const std::vector<double> &V, int p_threads) {
   
   long dim = (long)V.size();
   double val = 0;
   
#pragma omp parallel for reduction (+: val) num_threads (p_threads)
   for (long i = 0; i < dim; ++i) {
      val += V[i]*V[i];
   }
   
   return std::sqrt(val);
   
}
