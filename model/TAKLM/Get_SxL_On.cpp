#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SxL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   CRS Temp_SpL_On, Temp_SmL_On;
   Get_SpL_On(Temp_SpL_On, coeef*0.5);
   Get_SmL_On(Temp_SmL_On, coeef*0.5);
   Matrix_Matrix_Sum(Temp_SpL_On, Temp_SmL_On, M);
   Check_Symmetric_Matrix(M, zero_precision);
}
