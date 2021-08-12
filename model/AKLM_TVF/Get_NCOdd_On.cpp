#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_NCOdd_On(CRS &M, double coeef) {
   
   CRS COdd, COdd_D;
   Get_COdd_On(COdd, 1.0);
   Get_COdd_D_On(COdd_D, 1.0);
   Matrix_Matrix_Product(COdd_D, COdd, M);
   Matrix_Constant_Multiplication(M, coeef, 1);
   
}
