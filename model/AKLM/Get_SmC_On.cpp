#include "SML.hpp"
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_SmC_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SpC_On;
   Get_SpC_On(Temp_SpC_On, coeef);
   Matrix_Transpose(Temp_SpC_On, M);

}
