#include <iostream>
#include "SML.hpp"

void Change_Matrix_Basis(const CRS &M, const CCS &T, CRS &Out, int p_threads) {
   
   if (M.row_dim != M.col_dim || M.col_dim != T.row_dim) {
      std::cout << "Error in Change_Matrix_Basis" << std::endl;
      std::exit(0);
   }
   
   CCS Temp_Mat;
   std::vector<double> Temp_V1(std::max(T.col_dim, M.row_dim), 0.0);
   std::vector<double> Temp_V2(std::max(T.col_dim, M.row_dim), 0.0);
   
   double zero = 0.000000000000001;//pow(10,-15);
         
   //M*T --> Temp_Mat in CCS
   Temp_Mat.Col.push_back(0);
   for (int i = 0; i < T.col_dim; i++) {
      
#pragma omp parallel for num_threads (p_threads)
      for (long j = T.Col[i]; j < T.Col[i+1]; j++) {
         Temp_V1[T.Row[j]] = T.Val[j];
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int j = 0; j < M.row_dim; j++) {
         double val = 0;
         for (long k = M.Row[j]; k < M.Row[j+1]; k++) {
            val += Temp_V1[M.Col[k]]*M.Val[k];
         }
         Temp_V2[j] = val;
      }
      
      for (int j = 0; j < M.row_dim; j++) {
         if (std::abs(Temp_V2[j]) > zero) {
            Temp_Mat.Val.push_back(Temp_V2[j]);
            Temp_Mat.Row.push_back(j);
         }
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (long j = T.Col[i]; j < T.Col[i+1]; j++) {
         Temp_V1[T.Row[j]] = 0.0;
      }
      Temp_Mat.Col.push_back(Temp_Mat.Row.size());
   }
   
   Temp_Mat.row_dim = M.row_dim;
   Temp_Mat.col_dim = T.col_dim;
      
   //(T^†)*(M*T) --> Out in CRS
   Clear_Matrix(Out);
   Out.Row.push_back(0);
   for (int i = 0; i < T.col_dim; i++) {
      
#pragma omp parallel for num_threads (p_threads)
      for (long j = T.Col[i]; j < T.Col[i+1]; j++) {
         Temp_V1[T.Row[j]] = T.Val[j];
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int j = 0; j < Temp_Mat.col_dim; j++) {
         double val = 0;
         for (long k = Temp_Mat.Col[j]; k < Temp_Mat.Col[j+1]; k++) {
            val += Temp_V1[Temp_Mat.Row[k]]*Temp_Mat.Val[k];
         }
         Temp_V2[j] = val;
      }
      
      for (int j = 0; j < Temp_Mat.col_dim; j++) {
         if (std::abs(Temp_V2[j]) > zero) {
            Out.Val.push_back(Temp_V2[j]);
            Out.Col.push_back(j);
         }
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (long j = T.Col[i]; j < T.Col[i+1]; j++) {
         Temp_V1[T.Row[j]] = 0.0;
      }
      
      Out.Row.push_back(Out.Col.size());
      
   }
   
   Out.row_dim = T.col_dim;
   Out.col_dim = Temp_Mat.col_dim;

}
