#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SyC_2SyC_2_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_iSyC_2_On;
   Get_iSyC_2_On(Temp_iSyC_2_On, 1.0);
   Matrix_Matrix_Product(Temp_iSyC_2_On, Temp_iSyC_2_On, M);
   Matrix_Constant_Multiplication(M, -coeef, 1);
   Check_Symmetric_Matrix(M, zero_precision);
   
}
