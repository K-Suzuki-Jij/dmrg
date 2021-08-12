#include "SML.hpp"
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_SxC_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SpC_On, Temp_SmC_On;
   Get_SpC_On(Temp_SpC_On, coeef*0.5);
   Get_SmC_On(Temp_SmC_On, coeef*0.5);
   Matrix_Matrix_Sum(Temp_SpC_On, Temp_SmC_On, M);
   
   Check_Symmetric_Matrix(M, zero_precision);
}
