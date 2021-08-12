#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_NCDown_On(CRS &M, double coeef) {
   
   CRS CDown, CDown_D;
   Get_CDown_On(CDown, 1.0);
   Get_CDown_D_On(CDown_D, 1.0);
   Matrix_Matrix_Product(CDown_D, CDown, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
}
