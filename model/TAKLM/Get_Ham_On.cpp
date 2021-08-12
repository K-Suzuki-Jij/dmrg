#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Get_Ham_On(CRS &M) {
   
   Check_Parameters();
   
   CRS SxC_1, SxC_2, SxL, iSyC_1, iSyC_2, iSyL, SzC_1, SzC_2, SzL;
   
   Get_SxC_1_On (SxC_1 , 1.0);
   Get_SxC_2_On (SxC_2 , 1.0);
   Get_SxL_On   (SxL   , 1.0);
   Get_iSyC_1_On(iSyC_1, 1.0);
   Get_iSyC_2_On(iSyC_2, 1.0);
   Get_iSyL_On  (iSyL  , 1.0);
   Get_SzC_1_On (SzC_1 , 1.0);
   Get_SzC_2_On (SzC_2 , 1.0);
   Get_SzL_On   (SzL   , 1.0);
   
   CRS SxC_1SxL, SyC_1SyL, SzC_1SzL, SxC_2SxL, SyC_2SyL, SzC_2SzL;
   
   Matrix_Matrix_Product(SxC_1 , SxL , SxC_1SxL);
   Matrix_Matrix_Product(iSyC_1, iSyL, SyC_1SyL);
   Matrix_Matrix_Product(SzC_1 , SzL , SzC_1SzL);
   Matrix_Matrix_Product(SxC_2 , SxL , SxC_2SxL);
   Matrix_Matrix_Product(iSyC_2, iSyL, SyC_2SyL);
   Matrix_Matrix_Product(SzC_2 , SzL , SzC_2SzL);
   
   Matrix_Constant_Multiplication(SxC_1SxL, J_xy_1 , 1.0);
   Matrix_Constant_Multiplication(SyC_1SyL, -J_xy_1, 1.0);
   Matrix_Constant_Multiplication(SzC_1SzL, J_z_1  , 1.0);
   Matrix_Constant_Multiplication(SxC_2SxL, J_xy_2 , 1.0);
   Matrix_Constant_Multiplication(SyC_2SyL, -J_xy_2, 1.0);
   Matrix_Constant_Multiplication(SzC_2SzL, J_z_2  , 1.0);
   
   Matrix_Constant_Multiplication(SzC_1, -hc_z_1 , 1);
   Matrix_Constant_Multiplication(SzC_2, -hc_z_2 , 1);
   Matrix_Constant_Multiplication(SzL, -hl_z , 1);
   
   CRS SzLSzL;
   Get_SzLSzL_On(SzLSzL, -D_z);
   
   CRS NCUp_1, NCDown_1, NC_1, NCUp_1NCDown_1, NCUp_2, NCDown_2, NC_2, NCUp_2NCDown_2;
   
   Get_NCUp_1_On  (NCUp_1  , 1.0);
   Get_NCDown_1_On(NCDown_1, 1.0);
   Get_NC_1_On    (NC_1    , 1.0);
   Get_NCUp_2_On  (NCUp_2  , 1.0);
   Get_NCDown_2_On(NCDown_2, 1.0);
   Get_NC_2_On    (NC_2    , 1.0);
   Matrix_Matrix_Product(NCUp_1, NCDown_1, NCUp_1NCDown_1);
   Matrix_Matrix_Product(NCUp_2, NCDown_2, NCUp_2NCDown_2);

   Matrix_Constant_Multiplication(NCUp_1NCDown_1, U_1, 1);
   Matrix_Constant_Multiplication(NCUp_2NCDown_2, U_2, 1);
   Matrix_Constant_Multiplication(NC_1, -mu_1, 1);
   Matrix_Constant_Multiplication(NC_2, -mu_2, 1);

   CRS Temp1, Temp2;
   
   Matrix_Matrix_Sum(SxC_1SxL, SyC_1SyL      , Temp1);
   Matrix_Matrix_Sum(Temp1   , SzC_1SzL      , Temp2);
   Matrix_Matrix_Sum(Temp2   , SxC_2SxL      , Temp1);
   Matrix_Matrix_Sum(Temp1   , SyC_2SyL      , Temp2);
   Matrix_Matrix_Sum(Temp2   , SzC_2SzL      , Temp1);
   Matrix_Matrix_Sum(Temp1   , SzL           , Temp2);
   Matrix_Matrix_Sum(Temp2   , SzC_1         , Temp1);
   Matrix_Matrix_Sum(Temp1   , SzC_2         , Temp2);
   Matrix_Matrix_Sum(Temp2   , SzLSzL        , Temp1);
   Matrix_Matrix_Sum(Temp1   , NCUp_1NCDown_1, Temp2);
   Matrix_Matrix_Sum(Temp2   , NCUp_2NCDown_2, Temp1);
   Matrix_Matrix_Sum(Temp1   , NC_1          , Temp2);
   Matrix_Matrix_Sum(Temp2   , NC_2          , M    );

   Check_Symmetric_Matrix(M, zero_precision);
}
