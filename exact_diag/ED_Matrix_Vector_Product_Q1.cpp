#include <cmath>
#include <omp.h>
#include "ED.hpp"

void ED_Matrix_Vector_Product_Q1(CRS &M, std::vector<double> &Vec_In, std::vector<long> &Bases_In, std::vector<double> &Vec_Out, std::vector<long> &Bases_Out, int site, int p_threads) {
   
   if (Vec_In.size() != Bases_In.size() || Vec_Out.size() != Bases_Out.size()) {
      printf("Error in ED_Matrix_Vector_Product_Q1\n");
      exit(1);
   }
   
   int  dim_onsite     = M.row_dim;
   int  dim_target_in  = (int)Bases_In.size();
   int  dim_target_out = (int)Bases_Out.size();
   long site_constant  = std::pow(dim_onsite, site);
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < dim_target_out; i++) {
      long   basis_out       = Bases_Out[i];
      int    local_basis_out = ED_Find_Local_Basis(basis_out, site, dim_onsite);
      double temp_val        = 0;
      for (long j = M.Row[local_basis_out]; j < M.Row[local_basis_out + 1]; j++) {
         long a_basis_out = basis_out - (local_basis_out - M.Col[j])*site_constant;
         long inv_in      = Binary_Search(Bases_In, 0, dim_target_in, a_basis_out);
         if (inv_in >= 0) {
            temp_val += Vec_In[inv_in]*M.Val[j];
         }
      }
      Vec_Out[i] = temp_val;
   }
   
}
