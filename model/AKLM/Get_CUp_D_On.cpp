#include "SML.hpp"
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_CUp_D_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_CUp_On;
   Get_CUp_On(Temp_CUp_On, coeef);
   Matrix_Transpose(Temp_CUp_On, M);

}
