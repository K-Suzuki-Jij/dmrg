#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_SpC_On(CRS &M, double coeef) {
   
   CRS CUp_D, CDown;
   Get_CUp_D_On(CUp_D, 1.0);
   Get_CDown_On(CDown, 1.0);
   Matrix_Matrix_Product(CUp_D, CDown, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
}
