#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SpC_1_On(CRS &M, double coeef) {
   
   CRS CUp_1_D, CDown_1;
   Get_CUp_1_D_On(CUp_1_D, 1.0);
   Get_CDown_1_On(CDown_1, 1.0);
   Matrix_Matrix_Product(CUp_1_D, CDown_1, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
}
