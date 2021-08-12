#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_SzCSzC_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SzC_On;
   Get_SzC_On(Temp_SzC_On, 1.0);
   Matrix_Matrix_Product(Temp_SzC_On, Temp_SzC_On, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
   Check_Symmetric_Matrix(M, zero_precision);
   
}
