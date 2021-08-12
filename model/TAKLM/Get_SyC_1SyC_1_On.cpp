#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SyC_1SyC_1_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_iSyC_1_On;
   Get_iSyC_1_On(Temp_iSyC_1_On, 1.0);
   Matrix_Matrix_Product(Temp_iSyC_1_On, Temp_iSyC_1_On, M);
   Matrix_Constant_Multiplication(M, -coeef, 1);
   Check_Symmetric_Matrix(M, zero_precision, 1);
   
}
