#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SmC_1_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SpC_1_On;
   Get_SpC_1_On(Temp_SpC_1_On, coeef);
   Matrix_Transpose(Temp_SpC_1_On, M);

}
