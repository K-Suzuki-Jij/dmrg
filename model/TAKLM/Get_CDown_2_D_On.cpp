#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_CDown_2_D_On(CRS &M, double coeef) {
   
   CRS Temp_CDown_2_On;
   Get_CDown_2_On(Temp_CDown_2_On, coeef);
   Matrix_Transpose(Temp_CDown_2_On, M);

}
