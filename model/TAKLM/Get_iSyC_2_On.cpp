#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_iSyC_2_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SpC_2_On, Temp_SmC_2_On;
   Get_SpC_2_On(Temp_SpC_2_On, coeef*0.5);
   Get_SmC_2_On(Temp_SmC_2_On, -1.0*coeef*0.5);
   Matrix_Matrix_Sum(Temp_SpC_2_On, Temp_SmC_2_On, M);
   
}
