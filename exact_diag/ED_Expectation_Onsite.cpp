#include "ED.hpp"

void ED_Expectation_Onsite(CRS &M, std::vector<double> &Vec, std::vector<long> &Bases, std::vector<double> &Out, int tot_site, int p_threads) {
   
   if (Vec.size() != Bases.size() || (int)Out.size() != tot_site) {
      printf("Error in ED_Expectation_Onsite\n");
      exit(1);
   }
   
   std::vector<double> Temp_Vec(Bases.size(), 0.0);
   
   for (int site = 0; site < tot_site; site++) {
      ED_Matrix_Vector_Product_Q1(M, Vec, Bases, Temp_Vec, Bases, site, p_threads);
      Out[site] = Inner_Product(Vec, Temp_Vec, p_threads);
   }
   
}
