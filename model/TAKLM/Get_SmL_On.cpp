#include <cmath>
#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SmL_On(CRS &M, double coeef) {
   
  CRS Temp_SpL_On;
  Get_SpL_On(Temp_SpL_On, coeef);
  Matrix_Transpose(Temp_SpL_On, M);
   
}
