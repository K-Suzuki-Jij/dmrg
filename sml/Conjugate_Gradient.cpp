#include <iostream>
#include <iomanip>
#include <omp.h>
#include "SML.hpp"

void Conjugate_Gradient(const CRS &M, const std::vector<double> &Vec_In, std::vector<double> &Vec_Out, const Diag_Param &Diag_Param, int p_threads) {
   
   //Check input matrix
   if (M.row_dim != M.col_dim || M.row_dim <= 0 || M.col_dim <= 0) {
      std::cout << "Error in Conjugate_Gradient at 1" << std::endl;
      std::cout << "The input matrix is not a square one" << std::endl;
      std::cout << "row=" << M.row_dim << ", col=" << M.col_dim << std::endl;
      std::exit(0);
   }
   
   if (!(Vec_In.size() == Vec_Out.size() && (int)Vec_Out.size() == M.row_dim && (int)Vec_In.size() == M.row_dim)) {
      std::cout << "Error in Conjugate_Gradient at 2" << std::endl;
      std::cout << "Vec_Out_size=" << Vec_Out.size() << ", Vec_In_size" << Vec_In.size() << ", M_row_dim=" << M.row_dim << std::endl;
      std::exit(0);
   }
         
   double start = omp_get_wtime();
   double acc   = Diag_Param.cg_acc;
   int max_step = Diag_Param.cg_max_step;
   int row_dim  = M.row_dim;
   int step_num = 0;
   
   std::vector<double> RRR(row_dim, 0);
   std::vector<double> PPP(row_dim, 0);
   std::vector<double> YYY(row_dim, 0);
   std::vector<std::vector<double>> Temp_V;//Used only when Diag.Mat_Type = "Sym"
   
   if (Diag_Param.Mat_Type == "Sym") {
      Temp_V.resize(p_threads);
      for (int i = 0; i < p_threads; i++) {
         Temp_V[i].resize(row_dim);
      }
   }
   
   Normalize(Vec_Out, p_threads);
   
   Matrix_Vector_Product(M, Vec_Out, RRR, Temp_V, Diag_Param.Mat_Type, p_threads);
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < row_dim; i++) {
      RRR[i] = Vec_In[i] - RRR[i];
      PPP[i] = RRR[i];
   }
   
   for (int step = 1; step < max_step; step++) {
      
      Matrix_Vector_Product(M, PPP, YYY, Temp_V, Diag_Param.Mat_Type, p_threads);
      
      double inner_prod = Inner_Product(RRR, RRR, p_threads);
      double alpha      = inner_prod/Inner_Product(PPP, YYY, p_threads);
      
#pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < row_dim; i++) {
         Vec_Out[i] += alpha*PPP[i];
         RRR[i]     -= alpha*YYY[i];
      }
      
      double residual_error = Inner_Product(RRR, RRR, p_threads);

      std::cout << "\rCG_Step[" << step << "]=" << std::scientific << std::setprecision(1) << residual_error << std::string(5, ' ') << std::flush;
   
      if (residual_error < acc) {
         step_num = step;
         break;
      }
      
      double beta = residual_error/inner_prod;
      
#pragma omp parallel for num_threads (p_threads)
      for (int i = 0; i < row_dim; i++) {
         PPP[i] = RRR[i] + beta*PPP[i];
      }
   }
   
   std::cout << "\r";

   Output_Step_Number(step_num, omp_get_wtime() - start, "CG_Step.txt");
   
}
