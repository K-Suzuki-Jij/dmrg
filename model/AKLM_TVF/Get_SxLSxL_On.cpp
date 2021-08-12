#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_SxLSxL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SxL_On;
   Get_SxL_On(Temp_SxL_On, 1.0);
   Matrix_Matrix_Product(Temp_SxL_On, Temp_SxL_On, M);
   Matrix_Constant_Multiplication(M, coeef, 1);

   Check_Symmetric_Matrix(M, zero_precision, 1);
}

