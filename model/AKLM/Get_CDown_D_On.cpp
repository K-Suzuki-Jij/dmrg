#include "SML.hpp"
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_CDown_D_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_CDown_On;
   Get_CDown_On(Temp_CDown_On, coeef);
   Matrix_Transpose(Temp_CDown_On, M);

}
