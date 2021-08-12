#include "SML.hpp"
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_SzLSzL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SzL_On;
   Get_SzL_On(Temp_SzL_On, 1.0);
   Matrix_Matrix_Product(Temp_SzL_On, Temp_SzL_On, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
   Check_Symmetric_Matrix(M, zero_precision);
}
