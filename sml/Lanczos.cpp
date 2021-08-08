#include <iostream>
#include <iomanip>
#include <random>
#include <omp.h>
#include "SML.hpp"

void Lanczos(const CRS &M, std::vector<double> &Out_Vec, double &out_val, const Diag_Param &Diag_Param, int p_threads) {
   
   double start = omp_get_wtime();
   
   double   acc = Diag_Param.diag_acc;
   int min_step = Diag_Param.diag_min_step;
   int max_step = Diag_Param.diag_max_step;
   
   //Check input matrix
   if (M.row_dim != M.col_dim || M.row_dim <= 0 || M.col_dim <= 0) {
      std::cout << "Error in Lanczos at 1" << std::endl;
      std::cout << "The input matrix is not a square one" << std::endl;
      std::cout << "row=" << M.row_dim << ", col=" << M.col_dim << std::endl;
      std::exit(0);
   }
   
   if (M.row_dim == 1) {
      out_val = M.Val[0];
      Out_Vec.resize(M.row_dim);
      Out_Vec[0] = 1;
      return;
   }
   
   //Check step numbers
   if (max_step <= min_step) {
      std::cout << "Error in Lanczos at 2" << std::endl;
      std::cout << "min_step=" << min_step << ", max_step=" << max_step << std::endl;
      std::exit(0);
   }
   
   //Diagonalize by lapack if dim < 1000
   if (M.row_dim <= 1000) {
      Lapack_Dsyev(M, Out_Vec, out_val, 0);
      Output_Step_Number(0, omp_get_wtime() - start, "Lanczos_Step.txt");
      return;
   }
   
   int row_dim  = M.row_dim;
   int step_num = 0;
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
   if (Diag_Param.Lanczos_Initial_Guess == "Yes") {
      if (row_dim != (int)Out_Vec.size()) {
         std::cout << "Error in Lanczos at 3" << std::endl;
         std::exit(0);
      }
#pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < row_dim; i++) {
         Vector_0[i] = Out_Vec[i];
      }
   }
   else if (Diag_Param.Lanczos_Initial_Guess == "No") {
      mt.seed(seed);
      //Do not openmp here
      for (int i = 0; i < row_dim; i++) {
         Vector_0[i] = uniform_rand(mt);
      }
   }
   else {
      std::cout << "Error in Lanczos at 4" << std::endl;
      std::exit(0);
   }
   
   Normalize(Vector_0, p_threads);
   
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
      Matrix_Vector_Product(M, Vector_2, Vector_1, Temp_V, Diag_Param.Mat_Type, p_threads);
      Diag[step] = Inner_Product(Vector_1, Vector_2, p_threads);
      
      if (step >= min_step) {
         Lapack_Dstev(step + 1, Diag, Off_Diag, Temp_Eigen_Vec, Temp_Eigen_Val[step]);
         
         std::cout << "\rLanczos_Step[" << step << "]=" << std::scientific << std::setprecision(1);
         std::cout << std::abs(Temp_Eigen_Val[step] - Temp_Eigen_Val[step - 1]) << std::flush;
         
         if (std::abs(Temp_Eigen_Val[step] - Temp_Eigen_Val[step - 1]) < acc) {
            std::cout << "\r";
            step_num  = step;
            out_val = Temp_Eigen_Val[step];
            break;
         }
      }
      
#pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < row_dim; i++) {
         Vector_1[i] -= Diag[step]*Vector_2[i] + Off_Diag[step - 1]*Vector_0[i];
      }
      
      Copy_Vector(Vector_2, Vector_0, p_threads);
      
   }
   
   if (step_num == 0) {
      std::cout << "Error in Lanczos at 5" << std::endl;
      std::exit(0);
   }
   
   if (Diag_Param.Calc_Vec == "Yes") {
      //Restart Lanczos  //Restart Lanczos  //Restart Lanczos  //Restart Lanczos  //Restart Lanczos
      
      if ((int)Out_Vec.size() != row_dim) {
         Out_Vec.resize(row_dim);
      }
      
      //Set Initial Vector
      if (Diag_Param.Lanczos_Initial_Guess == "Yes") {
#pragma omp parallel for
         for (int i = 0; i < row_dim; i++) {
            Vector_0[i] = Out_Vec[i];
            Out_Vec[i]  = 0.0;
         }
      }
      else {
         mt.seed(seed);
         for (int i = 0; i < row_dim; i++) {
            Vector_0[i]  = uniform_rand(mt);
            Out_Vec[i] = 0.0;
         }
      }
      
      Normalize(Vector_0, p_threads);
      
#pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < row_dim; i++) {
         Out_Vec[i] += Temp_Eigen_Vec[0]*Vector_0[i];
      }
      
      Matrix_Vector_Product(M, Vector_0, Vector_1, Temp_V, Diag_Param.Mat_Type, p_threads);
      
#pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < row_dim; i++) {
         Vector_1[i] -= Diag[0]*Vector_0[i];
      }
      
      for (int step = 1; step <= step_num; step++) {
         
         Copy_Vector(Vector_1, Vector_2, p_threads);
         Normalize(Vector_2, p_threads);
         
#pragma omp parallel for num_threads (p_threads)
         for (int i = 0; i < row_dim; i++) {
            Out_Vec[i] += Temp_Eigen_Vec[step]*Vector_2[i];
         }
         
         Matrix_Vector_Product(M, Vector_2, Vector_1, Temp_V, Diag_Param.Mat_Type, p_threads);
         
#pragma omp parallel for num_threads (p_threads)
         for (int i = 0; i < row_dim; i++) {
            Vector_1[i] -= Diag[step]*Vector_2[i] + Off_Diag[step - 1]*Vector_0[i];
         }
         
         Copy_Vector(Vector_2, Vector_0, p_threads);
         
         std::cout << "\rLanczos_Vec_Step:" << step << "/" << step_num << std::string(5, ' ') << std::flush;
      }
      
      Normalize(Out_Vec, p_threads);
   }
   else if (Diag_Param.Calc_Vec != "No") {
      std::cout << "Error in Lanczos at 6" << std::endl;
      std::exit(0);
   }
   
   std::cout << "\r";
   
   Output_Step_Number(step_num, omp_get_wtime() - start, "Lanczos_Step.txt");
   
}

