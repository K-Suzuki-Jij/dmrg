//
//  Created by Kohei Suzuki on 2021/01/06.
//

#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_Ham_On(CRS &M) {
   
   dim_onsite = Find_Dim_Onsite();
   
   CRS ChemPot_mu;
   Get_NC_Tot_On(ChemPot_mu, ele_mu);
   
   CRS Mag_SzC_hz;
   Get_SzC_Tot_On(Mag_SzC_hz, ele_hz);
   
   CRS Mag_SzL_hz;
   Get_SzC_Tot_On(Mag_SzL_hz, lspin_hz);
   
   CRS Aniso_Dz;
   Get_SzLSzL_Tot_On(Aniso_Dz, lspin_Dz);
   
   CRS Coulomb_U, Temp1, Temp2;
   
   Temp2.row_dim = dim_onsite;
   Temp2.col_dim = dim_onsite;
   Coulomb_U.row_dim = dim_onsite;
   Coulomb_U.col_dim = dim_onsite;

   for (int i = 0; i <= Temp2.row_dim; i++) {
      Temp2.Row.push_back(0);
      Coulomb_U.Row.push_back(0);
   }
   
   for (int i = 0; i < num_ele_orbit; i++) {
      CRS Temp_NC_Up, Temp_NC_Down;
      Get_NC_Up_On  (Temp_NC_Up  , i, 1.0);
      Get_NC_Down_On(Temp_NC_Down, i, 1.0);
      Matrix_Matrix_Product(Temp_NC_Up, Temp_NC_Down, Temp1);
      Matrix_Constant_Multiplication(Temp1, ele_U, 1);
      Matrix_Matrix_Sum(Temp1, Temp2, Coulomb_U);
      Temp2 = Coulomb_U;
   }

   Clear_Matrix(Temp1);
   Clear_Matrix(Temp2);
   
   CRS Ele_Lspin_J_xy, Ele_Lspin_J_z;

   Temp1.row_dim = dim_onsite;
   Temp1.col_dim = dim_onsite;
   Temp2.row_dim = dim_onsite;
   Temp2.col_dim = dim_onsite;
   Ele_Lspin_J_xy.row_dim = dim_onsite;
   Ele_Lspin_J_xy.col_dim = dim_onsite;
   Ele_Lspin_J_z .row_dim = dim_onsite;
   Ele_Lspin_J_z .col_dim = dim_onsite;
   
   for (int i = 0; i <= Temp1.row_dim; i++) {
      Temp1.Row.push_back(0);
      Temp2.Row.push_back(0);
      Ele_Lspin_J_xy.Row.push_back(0);
      Ele_Lspin_J_z .Row.push_back(0);
   }
   
   for (int i = 0; i < num_ele_orbit; i++) {
      for (int j = 0; j < num_lspin_orbit; j++) {
         CRS Temp_SzC, Temp_SpC, Temp_SmC;
         CRS Temp_SzL, Temp_SpL, Temp_SmL;
         
         Get_SzC_On(Temp_SzC, i, 1.0);
         Get_SpC_On(Temp_SpC, i, 1.0);
         Get_SmC_On(Temp_SmC, i, 1.0);
         
         Get_SzL_On(Temp_SzL, j, 1.0);
         Get_SpL_On(Temp_SpL, j, 1.0);
         Get_SmL_On(Temp_SmL, j, 1.0);
         
         CRS Temp_SzCSzL, Temp_SpCSmL, Temp_SmCSpL;
         
         Matrix_Matrix_Product(Temp_SzC, Temp_SzL, Temp_SzCSzL);
         Matrix_Matrix_Product(Temp_SpC, Temp_SmL, Temp_SpCSmL);
         Matrix_Matrix_Product(Temp_SmC, Temp_SpL, Temp_SmCSpL);
         
         Matrix_Constant_Multiplication(Temp_SzCSzL, ele_lspin_Jz_orbit, 1);
         Matrix_Matrix_Sum(Temp_SzCSzL, Temp1, Ele_Lspin_J_z);
         Temp1 = Ele_Lspin_J_z;
         
         CRS Temp_SxCSxL_SyCSyL;
         Matrix_Matrix_Sum(Temp_SpCSmL, Temp_SmCSpL, Temp_SxCSxL_SyCSyL);
         Matrix_Constant_Multiplication(Temp_SxCSxL_SyCSyL, 0.5*ele_lspin_Jxy_orbit, 1);
         Matrix_Matrix_Sum(Temp_SxCSxL_SyCSyL, Temp2, Ele_Lspin_J_xy);
         Temp2 = Ele_Lspin_J_xy;
      }
   }

   Clear_Matrix(Temp1);
   Clear_Matrix(Temp2);
   
   CRS Lspin_Lspin_J_xy, Lspin_Lspin_J_z;
   
   Temp1.row_dim = dim_onsite;
   Temp1.col_dim = dim_onsite;
   Temp2.row_dim = dim_onsite;
   Temp2.col_dim = dim_onsite;
   Lspin_Lspin_J_xy.row_dim = dim_onsite;
   Lspin_Lspin_J_xy.col_dim = dim_onsite;
   Lspin_Lspin_J_z .row_dim = dim_onsite;
   Lspin_Lspin_J_z .col_dim = dim_onsite;
   
   for (int i = 0; i <= Temp1.row_dim; i++) {
      Temp1.Row.push_back(0);
      Temp2.Row.push_back(0);
      Lspin_Lspin_J_xy.Row.push_back(0);
      Lspin_Lspin_J_z .Row.push_back(0);
   }

   for (int i = 0; i < num_lspin_orbit - 1; i++) {
      CRS Temp_SzL_1, Temp_SpL_1, Temp_SmL_1;
      CRS Temp_SzL_2, Temp_SpL_2, Temp_SmL_2;

      Get_SzL_On(Temp_SzL_1, i  , 1.0);
      Get_SpL_On(Temp_SpL_1, i  , 1.0);
      Get_SmL_On(Temp_SmL_1, i  , 1.0);
      
      Get_SzL_On(Temp_SzL_2, i+1, 1.0);
      Get_SpL_On(Temp_SpL_2, i+1, 1.0);
      Get_SmL_On(Temp_SmL_2, i+1, 1.0);
      
      CRS Temp_SzL1SzL2, Temp_SpC1SmL2, Temp_SmC1SpL2;
      
      Matrix_Matrix_Product(Temp_SzL_1, Temp_SzL_2, Temp_SzL1SzL2);
      Matrix_Matrix_Product(Temp_SpL_1, Temp_SmL_2, Temp_SpC1SmL2);
      Matrix_Matrix_Product(Temp_SmL_1, Temp_SpL_2, Temp_SmC1SpL2);
      
      Matrix_Constant_Multiplication(Temp_SzL1SzL2, lspin_lspin_Jz_orbit, 1);
      
      Matrix_Matrix_Sum(Temp_SzL1SzL2, Temp1, Lspin_Lspin_J_z);
      Temp1 = Lspin_Lspin_J_z;
      
      CRS Temp_SxL1SxL2_SyL1SyL2;
      Matrix_Matrix_Sum(Temp_SpC1SmL2, Temp_SmC1SpL2, Temp_SxL1SxL2_SyL1SyL2);
      Matrix_Constant_Multiplication(Temp_SxL1SxL2_SyL1SyL2, 0.5*lspin_lspin_Jxy_orbit, 1);
      Matrix_Matrix_Sum(Temp_SxL1SxL2_SyL1SyL2, Temp2, Lspin_Lspin_J_xy);
      Temp2 = Lspin_Lspin_J_xy;
   }

   Clear_Matrix(Temp1);

   CRS Ele_Ele_V;
   
   Temp1.row_dim = dim_onsite;
   Temp1.col_dim = dim_onsite;
   Ele_Ele_V.col_dim = dim_onsite;
   Ele_Ele_V.row_dim = dim_onsite;
   
   for (int i = 0; i <= Temp1.row_dim; i++) {
      Temp1.Row.push_back(0);
      Ele_Ele_V.Row.push_back(0);
   }

   for (int i = 0; i < num_ele_orbit - 1; i++) {
      CRS Temp_NC_1, Temp_NC_2;
      Get_NC_On(Temp_NC_1, i  , 1);
      Get_NC_On(Temp_NC_2, i+1, 1);
      
      CRS Temp_NC1NC2;
      Matrix_Matrix_Product(Temp_NC_1, Temp_NC_2, Temp_NC1NC2);
      
      Matrix_Constant_Multiplication(Temp_NC1NC2, ele_ele_V_orbit, 1);
      
      Matrix_Matrix_Sum(Temp_NC1NC2, Temp1, Ele_Ele_V);
      Temp1 = Ele_Ele_V;
   }

   //Sum All of Avobe
   Matrix_Matrix_Sum(ChemPot_mu      , Mag_SzC_hz, Temp1);
   Matrix_Matrix_Sum(Mag_SzL_hz      , Temp1     , Temp2);
   Matrix_Matrix_Sum(Aniso_Dz        , Temp2     , Temp1);
   Matrix_Matrix_Sum(Coulomb_U       , Temp1     , Temp2);
   Matrix_Matrix_Sum(Ele_Lspin_J_xy  , Temp2     , Temp1);
   Matrix_Matrix_Sum(Ele_Lspin_J_z   , Temp1     , Temp2);
   Matrix_Matrix_Sum(Lspin_Lspin_J_xy, Temp2     , Temp1);
   Matrix_Matrix_Sum(Lspin_Lspin_J_z , Temp1     , Temp2);
   Matrix_Matrix_Sum(Ele_Ele_V       , Temp2     , M    );

}
