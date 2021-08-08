#include <iostream>
#include <vector>
#include <omp.h>

double Inner_Product(const std::vector<double> &V1, const std::vector<double> &V2, int p_threads) {
   
   if (V1.size() != V2.size()) {
      std::cout << "Error in Inner_Product\n" << std::endl;
      std::exit(0);
   }
   
   long   dim = (long)V1.size();
   double val = 0;
   
#pragma omp parallel for reduction (+: val) num_threads (p_threads)
   for (long i = 0; i < dim; ++i) {
      val += V1[i]*V2[i];
   }
   
   return val;
   
}
