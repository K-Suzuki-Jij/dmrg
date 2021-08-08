#include <iostream>
#include <vector>

void Copy_Vector(const std::vector<double> &V1, std::vector<double> &Copyed_V, int p_threads) {
   
   if (V1.size() != Copyed_V.size()) {
      std::cout << "Error in Copy_Vector\n" << std::endl;
      std::exit(0);
   }
   
   long dim = (long)V1.size();
   
#pragma omp parallel for num_threads (p_threads)
   for (long i = 0; i < dim; ++i) {
      Copyed_V[i] = V1[i];
   }
   
}
