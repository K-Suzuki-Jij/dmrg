#include <iostream>
#include <iomanip>
#include <random>
#include <omp.h>
#include "SML.hpp"

void Lanczos_Ex1(const CRS &M, const std::vector<double> &Eigen_Vec, double &out_val, const Diag_Param &Diag_Param, int p_threads) {
   
   double start = omp_get_wtime();
   
   double   acc = Diag_Param.diag_acc;
   int min_step = Diag_Param.diag_min_step;
   int max_step = Diag_Param.diag_max_step;
   
   //Check input matrix
   if (M.row_dim != M.col_dim || M.row_dim <= 0 || M.col_dim <= 0) {
      std::cout << "Error in Lanczos_Ex1 at 1" << std::endl;
      std::cout << "The input matrix is not a square one" << std::endl;
      std::cout << "row=" << M.row_dim << ", col=" << M.col_dim << std::endl;
      std::exit(0);
   }
   
   if (M.row_dim == 1) {
      out_val = M.Val[0];
      return;
   }
   
   //Check step numbers
   if (max_step <= min_step) {
      std::cout << "Error in Lanczos_Ex1 at 2" << std::endl;
      std::cout << "min_step=" << min_step << ", max_step=" << max_step << std::endl;
      std::exit(0);
   }
   
   if ((int)Eigen_Vec.size() != M.row_dim) {
      std::cout << "Error in Lanczos_Ex1 at 3" << std::endl;
      std::cout << "Eigen_Vec_size=" << Eigen_Vec.size() << ", row_dim=" << M.row_dim << std::endl;
      std::exit(0);
   }
   
   //Diagonalize by lapack if dim < 1000
   if (M.row_dim <= 1000) {
      std::vector<double> Temp_Out_Vec(M.row_dim, 0.0);
      Lapack_Dsyev(M, Temp_Out_Vec, out_val, 1);
      return;
   }
   
   int    row_dim  = M.row_dim;
   int    step_num = 0;
   double inn_pro  = 0.0;
   std::vector<double> Temp_Eigen_Val(max_step, 0);
   std::vector<double> Temp_Eigen_Vec;
   std::vector<double> Vector_0(row_dim, 0);
   std::vector<double> Vector_1(row_dim, 0);
   std::vector<double> Vector_2(row_dim, 0);
   std::vector<double> Diag(max_step, 0);
   std::vector<double> Off_Diag(max_step, 0);
   std::vector<std::vector<double>> Temp_V;//Used only when Diag.Mat_Type = "Sym"
   std::random_device rnd;
   std::uniform_real_distribution<> uniform_rand(-1, 1);
   std::mt19937 mt;
   int seed = rnd();
   
   if (Diag_Param.Mat_Type == "Sym") {
      Temp_V.resize(p_threads);
      for (int i = 0; i < p_threads; i++) {
         Temp_V[i].resize(row_dim);
#pragma omp parallel for num_threads (p_threads)
         for (int j = 0; j < row_dim; j++) {
            Temp_V[i][j] = 0.0;
         }
      }
   }
   
   //Set Initial Vector
   mt.seed(seed);
   for (int i = 0; i < row_dim; i++) {
      Vector_0[i] = uniform_rand(mt);
   }
   Normalize(Vector_0, p_threads);
   
   ///Orthoginalize  Orthoginalize  Orthoginalize
   inn_pro = Inner_Product(Vector_0, Eigen_Vec, p_threads);
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < row_dim; i++) {
      Vector_0[i] -= inn_pro*Eigen_Vec[i];
   }
   Normalize(Vector_0, p_threads);
   ///Orthoginalize  Orthoginalize  Orthoginalize
   
   Matrix_Vector_Product(M, Vector_0, Vector_1, Temp_V, Diag_Param.Mat_Type, p_threads);
   Diag[0] = Inner_Product(Vector_0, Vector_1, p_threads);
   Temp_Eigen_Val[0] = Diag[0];
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < row_dim; i++) {
      Vector_1[i] -= Diag[0]*Vector_0[i];
   }
   
   for (int step = 1; step < max_step; step++) {
      
      Copy_Vector(Vector_1, Vector_2, p_threads);
      Off_Diag[step - 1] = L2_Norm(Vector_2, p_threads);
      Normalize(Vector_2, p_threads);
      
      //////Orthoginalize  Orthoginalize  Orthoginalize//////
      inn_pro = Inner_Product(Vector_2, Eigen_Vec, p_threads);
#pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < row_dim; i++) {
         Vector_2[i] -= inn_pro*Eigen_Vec[i];
      }
      Normalize(Vector_2, p_threads);
      //////Orthoginalize  Orthoginalize  Orthoginalize//////
      
      Matrix_Vector_Product(M, Vector_2, Vector_1, Temp_V, Diag_Param.Mat_Type, p_threads);
      Diag[step] = Inner_Product(Vector_1, Vector_2, p_threads);
      
      if (step >= min_step) {
         Lapack_Dstev(step + 1, Diag, Off_Diag, Temp_Eigen_Vec, Temp_Eigen_Val[step]);
         
         std::cout << "\rLanczos_Ex1_Step[" << step << "]=" << std::scientific << std::setprecision(1);
         std::cout << std::abs(Temp_Eigen_Val[step] - Temp_Eigen_Val[step - 1]) << std::string(5, ' ') << std::flush;
         
         if (std::abs(Temp_Eigen_Val[step] - Temp_Eigen_Val[step - 1]) < acc) {
            step_num = step;
            out_val  = Temp_Eigen_Val[step];
            break;
         }
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < row_dim; i++) {
         Vector_1[i] -= Diag[step]*Vector_2[i] + Off_Diag[step - 1]*Vector_0[i];
      }
      
      Copy_Vector(Vector_2, Vector_0, p_threads);
      
   }

   std::cout << "\r";
   
   Output_Step_Number(step_num, omp_get_wtime() - start, "Lanczos_Ex1_Step.txt");
   
}
