#include <iostream>
#include <omp.h>
#include "SML.hpp"

void Matrix_Vector_Product(const CRS &M, const std::vector<double> &V, std::vector<double> &V_Out,  std::vector<std::vector<double>> &Temp_V, std::string Mat_Type, int p_threads) {
   
   if (Mat_Type == "Sym") {
      
      //Note that all the diagonal elements must be stored even if they are zero.
      if (M.row_dim != M.col_dim || (int)Temp_V.size() != p_threads || M.col_dim != (int)V.size()) {
         std::cout << "Error in Matrix_Vector_Product at 1" << std::endl;
         std::exit(0);
      }
      
      for (int thread_num = 0; thread_num < p_threads; thread_num++) {
         if((int)Temp_V[thread_num].size() != M.row_dim) {
            std::cout << "Error in Matrix_Vector_Product at 2" << std::endl;
            std::exit(0);
         }
      }
      
      if ((int)V_Out.size() != M.row_dim) {
         V_Out.resize(M.row_dim);
      }
      
      
   #pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < M.row_dim; ++i) {
         V_Out[i] = 0;
      }
      
   #pragma omp parallel for schedule (guided) num_threads (p_threads)
      for (long i = 0; i < M.row_dim; ++i) {
         int    thread_num = omp_get_thread_num();
         double temp_V     = V[i];
         double temp1      = M.Val[M.Row[i+1] - 1]*temp_V;
         for (long j = M.Row[i]; j < M.Row[i+1] - 1; ++j) {
            temp1 += M.Val[j]*V[M.Col[j]];
            Temp_V[thread_num][M.Col[j]] += M.Val[j]*temp_V;
         }
         Temp_V[thread_num][i] += temp1;
      }
      
   #pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < M.row_dim; ++i) {
         double temp = V_Out[i];
         for (int thread_num = 0; thread_num < p_threads; ++thread_num) {
            temp += Temp_V[thread_num][i];
            Temp_V[thread_num][i] = 0;
         }
         V_Out[i] = temp;
      }
      
      
   }
   else if (Mat_Type == "Non_Sym") {
      
      if (M.col_dim != (int)V.size()) {
         std::cout << "Error in Matrix_Vector_Product at 3" << std::endl;
         std::exit(0);
      }
      
      if (M.row_dim != (int)V_Out.size()) {
         V_Out.resize(M.row_dim);
      }

   #pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < M.row_dim; i++) {
         double temp = 0;
         for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
            temp += M.Val[j]*V[M.Col[j]];
         }
         V_Out[i] = temp;
      }
   }
   else {
      std::cout << "Error in Matrix_Vector_Product at 4" << std::endl;
      std::exit(0);
   }
   
}
