#include <iostream>
#include "SML.hpp"

extern "C" {
void dstev_(const char &JOBZ, const int &N, double *D, double *E, double **Z,  const int &LDZ, double* WORK, int& INFO);
};


void Lapack_Dstev(int dim, const std::vector<double> &Diag, const std::vector<double> &Off_Diag, std::vector<double> &Vec, double &val) {
   
   if ((int)Diag.size() < dim || (int)Off_Diag.size() < dim - 1) {
      std::cout << "Error in Lapack_Dstev" << std::endl;
      std::cout << "Diag_size=" << Diag.size() << ", Off_Diag_size=" << Off_Diag.size() << ", dim=" << dim << std::endl;
      std::exit(0);
   }
   
   if ((int)Vec.size() != dim) {
      Vec.resize(dim);
   }
   
   int info;
   double Lap_D[dim];
   double Lap_E[dim - 1];
   double Lap_Vec[dim][dim];
   double Lap_Work[2*dim];
   
   for (int i = 0; i < dim; i++) {
      Lap_D[i] = Diag[i];
   }
   
   for (int i = 0; i < dim - 1; i++) {
      Lap_E[i] = Off_Diag[i];
   }

   dstev_('V', dim, Lap_D, Lap_E, (double**)Lap_Vec, dim, Lap_Work, info);
   
   for (int i = 0; i < dim; i++) {
      Vec[i] = Lap_Vec[0][i];
   }
   
   val = Lap_D[0];
   
}
