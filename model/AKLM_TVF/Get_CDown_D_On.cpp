#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_CDown_D_On(CRS &M, double coeef) {
   
   CRS Temp;
   Get_CDown_On(Temp, coeef);
   Matrix_Transpose(Temp, M);
   
}
