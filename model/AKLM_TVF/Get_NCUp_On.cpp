#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_NCUp_On(CRS &M, double coeef) {
   
   CRS CUp, CUp_D;
   Get_CUp_On(CUp, 1.0);
   Get_CUp_D_On(CUp_D, 1.0);
   Matrix_Matrix_Product(CUp_D, CUp, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
}
