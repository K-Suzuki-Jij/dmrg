#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Get_Ham_On(CRS &M) {
   
   Check_Parameters();
   
   CRS SxC, SxL, iSyC, iSyL, SzC, SzL;
   
   Get_SxC_On (SxC, 1.0);
   Get_SxL_On (SxL, 1.0);
   Get_iSyC_On(iSyC, 1.0);
   Get_iSyL_On(iSyL, 1.0);
   Get_SzC_On (SzC, 1.0);
   Get_SzL_On (SzL, 1.0);
   
   CRS SxCSxL, SyCSyL, SzCSzL;
   
   Matrix_Matrix_Product(SxC, SxL, SxCSxL);
   Matrix_Matrix_Product(iSyC, iSyL, SyCSyL);
   Matrix_Matrix_Product(SzC, SzL, SzCSzL);
   
   Matrix_Constant_Multiplication(SxCSxL, J_xy , 1);
   Matrix_Constant_Multiplication(SyCSyL, -J_xy, 1);
   Matrix_Constant_Multiplication(SzCSzL, J_z , 1);
   
   Matrix_Constant_Multiplication(SzC, -hc_z , 1);
   Matrix_Constant_Multiplication(SzL, -hl_z , 1);
   
   CRS SzLSzL;
   Get_SzLSzL_On(SzLSzL, -D_z);
   
   CRS NCUp, NCDown, NC, NCUpNCDown;
   
   Get_NCUp_On  (NCUp  , 1.0);
   Get_NCDown_On(NCDown, 1.0);
   Get_NC_On    (NC    , 1.0);
   Matrix_Matrix_Product(NCUp, NCDown, NCUpNCDown);
   
   Matrix_Constant_Multiplication(NCUpNCDown, U, 1);
   Matrix_Constant_Multiplication(NC, -mu, 1);
   
   CRS Temp1, Temp2;
   
   Matrix_Matrix_Sum(SxCSxL, SyCSyL    , Temp1);
   Matrix_Matrix_Sum(Temp1 , SzCSzL    , Temp2);
   Matrix_Matrix_Sum(Temp2 , SzC       , Temp1);
   Check_Symmetric_Matrix(Temp1, zero_precision);

   Matrix_Matrix_Sum(Temp1 , SzL       , Temp2);
   Matrix_Matrix_Sum(Temp2 , SzLSzL    , Temp1);
   Matrix_Matrix_Sum(Temp1 , NCUpNCDown, Temp2);
   Matrix_Matrix_Sum(Temp2 , NC        , M    );
   
   Check_Symmetric_Matrix(M, zero_precision);
}
