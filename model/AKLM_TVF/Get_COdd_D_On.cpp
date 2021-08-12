#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_COdd_D_On(CRS &M, double coeef) {
   
   CRS Temp;
   Get_COdd_On(Temp, coeef);
   Matrix_Transpose(Temp, M);
   
}
