#include "ED.hpp"

void ED_Expectation_Intersite_Q1(CRS &M, std::vector<double> &Vec_In, std::vector<long> &Bases_In, std::vector<long> &Bases_Out, std::vector<double> &Out, int site_ref, int tot_site, int p_threads) {
   
   if (Vec_In.size() != Bases_In.size() || (int)Out.size() != tot_site - site_ref) {
      printf("Error in ED_Expectation_Intersite_Q1\n");
      exit(1);
   }
   
   std::vector<double> Temp_Vec_Refe(Bases_Out.size(), 0.0);
   std::vector<double> Temp_Vec_Dist(Bases_Out.size(), 0.0);
   
   ED_Matrix_Vector_Product_Q1(M, Vec_In, Bases_In, Temp_Vec_Refe, Bases_Out, site_ref, p_threads);

   for (int site = site_ref; site < tot_site; site++) {
      ED_Matrix_Vector_Product_Q1(M, Vec_In, Bases_In, Temp_Vec_Dist, Bases_Out, site, p_threads);
      Out[site - site_ref] = Inner_Product(Temp_Vec_Refe, Temp_Vec_Dist, p_threads);
   }
   
}

void ED_Expectation_Intersite_Q1(CRS &M, std::vector<double> &Vec_In, std::vector<long> &Bases_In, std::vector<long> &Bases_Out_1, std::vector<long> &Bases_Out_2, std::vector<double> &Out, int site_ref, int tot_site, int p_threads) {
   
   if (Vec_In.size() != Bases_In.size() || (int)Out.size() != tot_site - site_ref) {
      printf("Error in ED_Expectation_Intersite_Q1\n");
      exit(1);
   }
   
   std::vector<double> Temp_Vec_Refe(Bases_Out_1.size(), 0.0);
   std::vector<double> Temp_Vec_Dist(Bases_Out_1.size(), 0.0);
   
   ED_Matrix_Vector_Product_Q1(M, Vec_In, Bases_In, Temp_Vec_Refe, Bases_Out_1, site_ref, p_threads);

   for (int site = site_ref; site < tot_site; site++) {
      ED_Matrix_Vector_Product_Q1(M, Vec_In, Bases_In, Temp_Vec_Dist, Bases_Out_1, site, p_threads);
      Out[site - site_ref] = Inner_Product(Temp_Vec_Refe, Temp_Vec_Dist, p_threads);
   }
   
   Temp_Vec_Refe.resize(Bases_Out_2.size());
   Temp_Vec_Dist.resize(Bases_Out_2.size());

   ED_Matrix_Vector_Product_Q1(M, Vec_In, Bases_In, Temp_Vec_Refe, Bases_Out_2, site_ref, p_threads);

   for (int site = site_ref; site < tot_site; site++) {
      ED_Matrix_Vector_Product_Q1(M, Vec_In, Bases_In, Temp_Vec_Dist, Bases_Out_2, site, p_threads);
      Out[site - site_ref] += Inner_Product(Temp_Vec_Refe, Temp_Vec_Dist, p_threads);
   }
   
}
