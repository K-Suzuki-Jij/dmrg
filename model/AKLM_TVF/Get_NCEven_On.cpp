#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_NCEven_On(CRS &M, double coeef) {
   
   CRS CEven, CEven_D;
   Get_CEven_On(CEven, 1.0);
   Get_CEven_D_On(CEven_D, 1.0);
   Matrix_Matrix_Product(CEven_D, CEven, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
}
