#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SC_2SL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   
   CRS SzC_2, SzL, SpL, SpC_2, SmL, SmC_2;
   
   Get_SzC_2_On(SzC_2, 1.0);
   Get_SzL_On(SzL, 1.0);
   Get_SpC_2_On(SpC_2, 1.0);
   Get_SpL_On(SpL, 1.0);
   Get_SmC_2_On(SmC_2, 1.0);
   Get_SmL_On(SmL, 1.0);
   
   CRS SzC_2SzL, SpC_2SmL, SmC_2SpL;
   
   Matrix_Matrix_Product(SzC_2, SzL, SzC_2SzL);
   Matrix_Matrix_Product(SpC_2, SmL, SpC_2SmL);
   Matrix_Matrix_Product(SmC_2, SpL, SmC_2SpL);
   
   Matrix_Constant_Multiplication(SpC_2SmL, 0.5 , 1);
   Matrix_Constant_Multiplication(SmC_2SpL, 0.5 , 1);

   CRS Temp1;
   
   Matrix_Matrix_Sum(SmC_2SpL, SpC_2SmL , Temp1);
   Matrix_Matrix_Sum(Temp1   , SzC_2SzL , M    );
   
   Matrix_Constant_Multiplication(M, coeef, 1);
 
   Check_Symmetric_Matrix(M, zero_precision);
}
