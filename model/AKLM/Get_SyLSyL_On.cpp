#include "SML.hpp"
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_SyLSyL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_iSyL_On;
   Get_iSyL_On(Temp_iSyL_On, 1.0);
   Matrix_Matrix_Product(Temp_iSyL_On, Temp_iSyL_On, M);
   Matrix_Constant_Multiplication(M, -coeef, 1);

   Check_Symmetric_Matrix(M, zero_precision);
}
