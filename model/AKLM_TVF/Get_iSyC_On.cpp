#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_iSyC_On(CRS &M, double coeef) {
   
   CRS Temp_SpC_On, Temp_SmC_On;
   Get_SpC_On(Temp_SpC_On, coeef*0.5);
   Get_SmC_On(Temp_SmC_On, -1.0*coeef*0.5);
   Matrix_Matrix_Sum(Temp_SpC_On, Temp_SmC_On, M);
   
}
