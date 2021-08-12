#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_SyCSyC_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_iSyC_On;
   Get_iSyC_On(Temp_iSyC_On, 1.0);
   Matrix_Matrix_Product(Temp_iSyC_On, Temp_iSyC_On, M);
   Matrix_Constant_Multiplication(M, -coeef, 1);

   Check_Symmetric_Matrix(M, zero_precision, 1);
}
