#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SzC_2SzC_2_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SzC_2_On;
   Get_SzC_2_On(Temp_SzC_2_On, 1.0);
   Matrix_Matrix_Product(Temp_SzC_2_On, Temp_SzC_2_On, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   Check_Symmetric_Matrix(M, zero_precision);
   
}
