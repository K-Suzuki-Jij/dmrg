#include "SML.hpp"
#include "Model_1D_AKLM_TVF.hpp"

void Model_1D_AKLM_TVF::Get_SxC_On(CRS &M, double coeef) {
   
   CRS NCEven, NCOdd;
   Get_NCEven_On(NCEven, 1.0);
   Get_NCOdd_On(NCOdd, -1.0);
   Matrix_Matrix_Sum(NCEven, NCOdd, M);
   Matrix_Constant_Multiplication(M, coeef*0.5, 1);
   Check_Symmetric_Matrix(M, zero_precision);
   
}
