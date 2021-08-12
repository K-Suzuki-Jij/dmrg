#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_CDown_1_D_On(CRS &M, double coeef) {
   
   CRS Temp_CDown_1_On;
   Get_CDown_1_On(Temp_CDown_1_On, coeef);
   Matrix_Transpose(Temp_CDown_1_On, M);

}
