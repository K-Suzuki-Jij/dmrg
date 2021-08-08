#include <iostream>
#include "SML.hpp"

extern "C" {
void dsyev_(const char &JOBZ, const char &UPLO, const int &N, double **A, const int &LDA, double *W, double *WORK, const int &LWORK, int &INFO);
};


void Lapack_Dsyev(const CRS &M, std::vector<double> &Vec, double &val, int state) {
   
   if (M.row_dim != M.col_dim || M.row_dim == 0 || M.col_dim == 0) {
      std::cout << "Error in Lapack_Dsyev" << std::endl;
      std::cout << "The input matrix is not a square one" << std::endl;
      std::cout << "row=" << M.row_dim << ", col=" << M.col_dim << std::endl;
      std::exit(0);
   }
   
   if (state < 0 || state >= M.row_dim) {
      std::cout << "Error in Lapack_Dsyev\n" << std::endl;
      std::exit(0);
   }
   
   if ((int)Vec.size() != M.row_dim) {
      Vec.resize(M.row_dim);
   }
   
   double M_Array[M.row_dim][M.row_dim];
   
   //Initialize Arrays
   for (int i = 0; i < M.row_dim; i++) {
      for (int j = 0; j < M.row_dim; j++) {
         M_Array[i][j] = 0.0;
      }
   }
   
   //Input Matrix
   for (int i = 0; i < M.row_dim; i++) {
      for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
         M_Array[i][M.Col[j]] = M.Val[j];
         M_Array[M.Col[j]][i] = M.Val[j];
      }
   }
   
   //Lapack dsyev
   int info;
   double Val_Array[M.row_dim];
   double Work[3*M.row_dim];
   
   dsyev_('V', 'L', M.row_dim, (double**)M_Array, M.row_dim, Val_Array, Work, 3*M.row_dim, info);

   for (int i = 0; i < M.row_dim; i++) {
      Vec[i] = M_Array[state][i];
   }
   
   val = Val_Array[state];
   
}

void Lapack_Dsyev(const std::vector<std::vector<double>> &M, std::vector<std::vector<double>> &Vec, std::vector<double> &Val) {
   
   int row_dim = (int)M.size();
   
   for (int i = 0; i < row_dim; i++) {
      if (row_dim != (int)M[i].size()) {
         std::cout << "Error in Lapack_Dsyev" << std::endl;
         std::exit(0);
      }
   }
   
   
   if (row_dim != (int)Vec.size()) {
      Vec.resize(row_dim);
   }
   
   for (int i = 0; i < row_dim; i++) {
      if (row_dim != (int)Vec[i].size()) {
         Vec[i].resize(row_dim);
      }
   }
   
   if (row_dim != (int)Val.size()) {
      Val.resize(row_dim);
   }
   
   
   double M_Array[row_dim][row_dim];
   
   //Initialize Arrays
   for (int i = 0; i < row_dim; i++) {
      for (int j = 0; j < row_dim; j++) {
         M_Array[i][j] = 0.0;
      }
   }
   
   //Input Matrix
   for (int i = 0; i < row_dim; i++) {
      for (int j = i; j < row_dim; j++) {
         M_Array[i][j] = M[i][j];
         M_Array[j][i] = M[i][j];
      }
   }
   
   //Lapack dsyev
   int info;
   double Val_Array[row_dim];
   double Work[3*row_dim];
   
   dsyev_('V', 'L', row_dim, (double**)M_Array, row_dim, Val_Array, Work, 3*row_dim, info);

   for (int i = 0; i < row_dim; i++) {
      for (int j = 0; j < row_dim; j++) {
         Vec[i][j] = M_Array[i][j];
      }
      Val[i] = Val_Array[i];
   }
      
}

