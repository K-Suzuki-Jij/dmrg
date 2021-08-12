#include <cmath>
#include "SML.hpp"
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_NC_On(CRS &M, double coeef) {
   
   CRS Temp_NC_1, Temp_NC_2;
   Get_NC_1_On(Temp_NC_1, coeef);
   Get_NC_2_On(Temp_NC_2, coeef);
   Matrix_Matrix_Sum(Temp_NC_1, Temp_NC_2, M);
   Check_Symmetric_Matrix(M, zero_precision, 1);
   
}
