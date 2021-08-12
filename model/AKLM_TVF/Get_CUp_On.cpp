#include <cmath>
#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_CUp_On(CRS &M, double coeef) {
   
   CRS CEven, COdd;
   Get_CEven_On(CEven, 1.0);
   Get_COdd_On(COdd, 1.0);
   Matrix_Matrix_Sum(CEven, COdd, M);
   Matrix_Constant_Multiplication(M, coeef*(1.0/std::sqrt(2)), 1);
   
}
