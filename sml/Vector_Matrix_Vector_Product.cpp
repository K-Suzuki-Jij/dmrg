#include <iostream>
#include "SML.hpp"

double Vector_Matrix_Vector_Product(std::vector<double> &V1, CRS &M, std::vector<double> &V2, std::string Mat_Type, int p_threads) {
   
   if ((int)V1.size() != M.row_dim || (int)V2.size() != M.col_dim) {
      printf("Error in Vector_Matrix_Vector_Product\n");
      printf("v1_size=%lu, M_row_dim=%d, v2_size=%lu, M_col_dim=%d\n", V1.size(), M.row_dim, V2.size(), M.col_dim);
      exit(1);
   }
   
   double temp = 0;
   
   if (Mat_Type == "Non_Sym") {
#pragma omp parallel for reduction (+:temp) num_threads (p_threads)
      for (int i = 0; i < M.row_dim; i++) {
         double v1_val = V1[i];
         for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
            temp += v1_val*M.Val[j]*V2[M.Col[j]];
         }
      }
   }
   else if (Mat_Type == "Sym") {
      std::vector<std::vector<double>> Temp_V(p_threads, std::vector<double>(M.row_dim, 0.0));
      std::vector<double> Out_V(M.row_dim, 0.0);
      Matrix_Vector_Product(M, V2, Out_V, Temp_V, Mat_Type, p_threads);
      temp = Inner_Product(V1, V2, p_threads);
   }
   else {
      std::cout << "Error in Vector_Matrix_Vector_Product" << std::endl;
      exit(1);
   }
   
   return temp;
   
}
