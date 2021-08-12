#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_SC_1SL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   
   CRS SzC_1, SzL, SpL, SpC_1, SmL, SmC_1;
   
   Get_SzC_1_On(SzC_1, 1.0);
   Get_SzL_On(SzL, 1.0);
   Get_SpC_1_On(SpC_1, 1.0);
   Get_SpL_On(SpL, 1.0);
   Get_SmC_1_On(SmC_1, 1.0);
   Get_SmL_On(SmL, 1.0);
   
   CRS SzC_1SzL, SpC_1SmL, SmC_1SpL;
   
   Matrix_Matrix_Product(SzC_1, SzL, SzC_1SzL);
   Matrix_Matrix_Product(SpC_1, SmL, SpC_1SmL);
   Matrix_Matrix_Product(SmC_1, SpL, SmC_1SpL);
   
   Matrix_Constant_Multiplication(SpC_1SmL, 0.5 , 1);
   Matrix_Constant_Multiplication(SmC_1SpL, 0.5 , 1);

   CRS Temp1;
   
   Matrix_Matrix_Sum(SmC_1SpL, SpC_1SmL , Temp1);
   Matrix_Matrix_Sum(Temp1   , SzC_1SzL , M    );
   
   Matrix_Constant_Multiplication(M, coeef, 1);
 
   Check_Symmetric_Matrix(M, zero_precision);
}
