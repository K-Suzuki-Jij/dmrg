#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_SCSL_On(CRS &M, double coeef) {
   
   Check_Parameters();
   
   CRS SzC, SzL, SpL, SpC, SmL, SmC;
   
   Get_SzC_On(SzC, 1.0);
   Get_SzL_On(SzL, 1.0);
   Get_SpL_On(SpL, 1.0);
   Get_SpC_On(SpC, 1.0);
   Get_SmL_On(SmL, 1.0);
   Get_SmC_On(SmC, 1.0);
   
   CRS SzCSzL, SpCSmL, SmCSpL;
   
   Matrix_Matrix_Product(SzC, SzL, SzCSzL);
   Matrix_Matrix_Product(SpC, SmL, SpCSmL);
   Matrix_Matrix_Product(SmC, SpL, SmCSpL);
   
   Matrix_Constant_Multiplication(SpCSmL, 0.5 , 1);
   Matrix_Constant_Multiplication(SmCSpL, 0.5 , 1);

   CRS Temp1;
   
   Matrix_Matrix_Sum(SmCSpL, SpCSmL , Temp1);
   Matrix_Matrix_Sum(Temp1 , SzCSzL , M    );
   
   Matrix_Constant_Multiplication(M, coeef, 1);
 
   Check_Symmetric_Matrix(M, zero_precision);
}
