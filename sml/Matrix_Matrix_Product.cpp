#include <iostream>
#include <cmath>
#include "SML.hpp"

void Matrix_Matrix_Product(const CRS &M1, const CRS &M2, CRS &Out) {
   
   if (M1.col_dim != M2.row_dim) {
      std::cout << "Error in Matrix_Matrix_Product" << std::endl;
      std::cout << "Matrix product cannot be defined" << std::endl;
      std::cout << "M1.col_dim=" << M1.col_dim << ", M2.row_dim=" << M2.row_dim << std::endl;
      std::exit(0);
   }
   
   Clear_Matrix(Out);
   
   std::vector<double> Temp_V1(M1.col_dim, 0);
   std::vector<double> Temp_V2(M2.col_dim, 0);
   double zero = 0.000000000000001;//pow(10,-15);

   Out.Row.push_back(0);

   for (int i = 0; i < M1.row_dim; i++) {
      
      for (long j = M1.Row[i]; j < M1.Row[i+1]; j++) {
         Temp_V1[M1.Col[j]] = M1.Val[j];
      }
      
      for (int j = 0; j < M1.col_dim; j++) {
         for (long k = M2.Row[j]; k < M2.Row[j+1]; k++) {
            Temp_V2[M2.Col[k]] += Temp_V1[j]*M2.Val[k];
         }
      }
      
      for (int j = 0; j < M2.col_dim; j++) {
         if (std::abs(Temp_V2[j]) > zero) {
            Out.Val.push_back(Temp_V2[j]);
            Out.Col.push_back(j);
         }
      }
      
      for (long j = M1.Row[i]; j < M1.Row[i+1]; j++) {
         Temp_V1[M1.Col[j]] = 0;
      }
      
      for (int j = 0; j < M2.col_dim; j++) {
         Temp_V2[j] = 0;
      }
      Out.Row.push_back(Out.Col.size());
   }
   
   
   Out.row_dim = M1.row_dim;
   Out.col_dim = M2.col_dim;
   
}































