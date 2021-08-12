#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SpC_2_On(CRS &M, double coeef) {
   
   CRS CUp_2_D, CDown_2;
   Get_CUp_2_D_On(CUp_2_D, 1.0);
   Get_CDown_2_On(CDown_2, 1.0);
   Matrix_Matrix_Product(CUp_2_D, CDown_2, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
}
