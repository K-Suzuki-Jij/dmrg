#include <iostream>
#include <omp.h>
#include "SML.hpp"

void Inverse_Iteration(CRS &M, std::vector<double> &Eigen_Vec, double eigen_val, const Diag_Param &Diag_Param, int p_threads) {
   
   //Check input matrix
   if (M.row_dim != M.col_dim || M.row_dim <= 0 || M.col_dim <= 0) {
      std::cout << "Error in Inverse_Iteration" << std::endl;
      std::cout << "The input matrix is not a square one" << std::endl;
      std::cout << "row_dim=" << M.row_dim << ", col_dim=" << M.col_dim << std::endl;
      std::exit(0);
   }
   
   if ((int)Eigen_Vec.size() != M.row_dim) {
      std::cout << "Error in Inverse_Iteration" << std::endl;
      std::cout << "Eigen_Vec_size=" << Eigen_Vec.size() << ", M_row_dim=" << M.row_dim << std::endl;
      std::exit(0);
   }
   
   double start = omp_get_wtime();
   
   double diag_add = Diag_Param.ii_diag_add;
   double ii_acc   = Diag_Param.ii_acc;
   int ii_max_step = Diag_Param.ii_max_step;
   int step_num    = 0;
   
   std::vector<double> Improved_Eigen_Vec = Eigen_Vec;
   
   Diag_Add_Matrix(M, diag_add - eigen_val, p_threads);
   
   for (int step = 0; step < ii_max_step; step++) {
            
      double residual_error = Residual_Error_Eigenpair(M, Eigen_Vec, diag_add, Diag_Param.Mat_Type, p_threads);

      if (residual_error < ii_acc) {
         step_num = step;
         break;
      }
      
      if (step == ii_max_step) {
         std::cout << "Inverse_Iteration Not Converge, error=" << residual_error << std::endl;
         std::cout << "Continue..." << std::endl;
         break;
      }
      
      Conjugate_Gradient(M, Eigen_Vec, Improved_Eigen_Vec, Diag_Param, p_threads);
      
      Normalize(Improved_Eigen_Vec, p_threads);
      
      Copy_Vector(Improved_Eigen_Vec, Eigen_Vec, p_threads);
      
   }
   
   Diag_Add_Matrix(M, -(diag_add - eigen_val), p_threads);
   
   Output_Step_Number(step_num, omp_get_wtime() - start, "II_Step.txt");
   
}
