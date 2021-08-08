#include <iostream>
#include <cmath>
#include "SML.hpp"

double Residual_Error_Eigenpair(const CRS &M, const std::vector<double> &Eigen_Vec, double eigen_val, std::string Mat_Type, int p_threads) {
   
   double norm = 0;
   
   if (Mat_Type == "Non_Sym") {
#pragma omp parallel for reduction (+:norm) num_threads (p_threads)
      for (int i = 0; i < M.row_dim; i++) {
         double temp = 0;
         for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
            temp += M.Val[j]*Eigen_Vec[M.Col[j]];
         }
         norm += std::abs(temp - eigen_val*Eigen_Vec[i]);
      }
   }
   else if (Mat_Type == "Sym") {
      std::vector<double> Temp_V(M.row_dim, 0);
      for (int i = 0; i < M.row_dim; i++) {
         double temp1 = M.Val[M.Row[i + 1] - 1]*Eigen_Vec[i];
         double temp2 = Eigen_Vec[i];
         for (long j = M.Row[i]; j < M.Row[i+1] - 1; j++)  {
            temp1 += M.Val[j]*Eigen_Vec[M.Col[j]];
            Temp_V[M.Col[j]] += M.Val[j]*temp2;
         }
         Temp_V[i] += temp1;
      }
#pragma omp parallel for reduction (+:norm) num_threads (p_threads)
      for (int i = 0; i < M.row_dim; i++) {
         norm += std::abs(eigen_val*Eigen_Vec[i] - Temp_V[i]);
      }
   }
   else {
      std::cout << "Error in Residual_Error_Eigenpair" << std::endl;
      std::exit(0);
   }
   
   return norm;
   
}
