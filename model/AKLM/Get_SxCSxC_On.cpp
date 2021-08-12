#include "SML.hpp"
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_SxCSxC_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SxC_On;
   Get_SxC_On(Temp_SxC_On, 1.0);
   Matrix_Matrix_Product(Temp_SxC_On, Temp_SxC_On, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
   Check_Symmetric_Matrix(M, zero_precision);
}
