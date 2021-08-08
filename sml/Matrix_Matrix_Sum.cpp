#include <iostream>
#include <cmath>
#include "SML.hpp"

void Matrix_Matrix_Sum(const CRS &M1, const CRS &M2, CRS &Out) {
   
   if (M1.row_dim != M2.row_dim || M1.col_dim != M2.col_dim) {
      std::cout << "Error in Matrix_Matrix_Sum" << std::endl;
      std::cout << "Matrix types do not match each other" << std::endl;
      
      printf("M1.row_dim=%d\n", M1.row_dim);
      printf("M2.row_dim=%d\n", M2.row_dim);
      printf("M1.col_dim=%d\n", M1.col_dim);
      printf("M2.col_dim=%d\n", M2.col_dim);
      exit(1);
   }
   
   double zero = 0.000000000000001;//pow(10,-15);
   
   std::vector<std::vector<double>> MM1(M1.row_dim, std::vector<double>(M1.col_dim, 0.0));
   std::vector<std::vector<double>> MM2(M1.row_dim, std::vector<double>(M1.col_dim, 0.0));
   
   Clear_Matrix(Out);
   Out.Row.push_back(0);
   for (int i = 0; i < M1.col_dim; i++) {
      for (long j = M1.Row[i]; j < M1.Row[i+1]; j++) {
         MM1[i][M1.Col[j]] = M1.Val[j];
      }
   }
   for (int i = 0; i < M2.col_dim; i++) {
      for (long j = M2.Row[i]; j < M2.Row[i+1]; j++) {
         MM2[i][M2.Col[j]] = M2.Val[j];
      }
   }

   for (int i = 0; i < M1.row_dim; i++) {
      for (int j = 0; j < M1.col_dim; j++) {
         if (std::abs(MM1[i][j] + MM2[i][j]) > zero) {
            Out.Val.push_back(MM1[i][j] + MM2[i][j]);
            Out.Col.push_back(j);
         }
      }
      Out.Row.push_back(Out.Col.size());
   }
   
   Out.row_dim = M1.row_dim;
   Out.col_dim = M1.col_dim;
   
}
