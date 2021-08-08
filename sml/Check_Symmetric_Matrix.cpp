#include <iostream>
#include "SML.hpp"

void Check_Symmetric_Matrix(const CRS &M, double zero_precision) {
   
   for (int i = 0; i < M.row_dim; i++) {
      for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
         long inv = Binary_Search(M.Col, M.Row[M.Col[j]], M.Row[M.Col[j] + 1], i);
         if (inv == -1) {
            std::cout << "The input matrix is not symmetric." << std::endl;
            std::cout << "Corresponding element does not exist." << std::endl;
            std::cout << "row=" << i << ", col=" << M.Col[j] << ", val=" << M.Val[j] << std::endl;
            std::exit(0);
         }
         if (std::abs(M.Val[j] - M.Val[inv]) > zero_precision) {
            std::cout << "The input matrix is not symmetric." << std::endl;
            std::cout << "M[" << i << "][" << M.Col[j] << "]=" << M.Val[j] << ", " << M.Val[inv] << "=M[" << M.Col[j] << "][" << i << "]" << std::endl;
            std::exit(0);
         }
      }
   }
}
