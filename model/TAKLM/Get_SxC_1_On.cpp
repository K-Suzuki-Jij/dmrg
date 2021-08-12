#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SxC_1_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SpC_1_On, Temp_SmC_1_On;
   Get_SpC_1_On(Temp_SpC_1_On, coeef*0.5);
   Get_SmC_1_On(Temp_SmC_1_On, coeef*0.5);
   Matrix_Matrix_Sum(Temp_SpC_1_On, Temp_SmC_1_On, M);
   Check_Symmetric_Matrix(M, zero_precision);
   
}
