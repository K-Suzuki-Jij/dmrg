//
//  Created by Kohei Suzuki on 2021/01/06.
//

#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_SCSL_On(CRS &M, int target_ele_orbit, int target_lspin_orbit, double coeef) {
   
   CRS Temp_SzC, Temp_SpC, Temp_SmC;
   CRS Temp_SzL, Temp_SpL, Temp_SmL;
   
   Get_SzC_On(Temp_SzC, target_ele_orbit  , 1.0);
   Get_SpC_On(Temp_SpC, target_ele_orbit  , 1.0);
   Get_SmC_On(Temp_SmC, target_ele_orbit  , 1.0);
   
   Get_SzL_On(Temp_SzL, target_lspin_orbit, 1.0);
   Get_SpL_On(Temp_SpL, target_lspin_orbit, 1.0);
   Get_SmL_On(Temp_SmL, target_lspin_orbit, 1.0);
   
   CRS Temp_SzCSzL, Temp_SpCSmL, Temp_SmCSpL;
   
   Matrix_Matrix_Product(Temp_SzC, Temp_SzL, Temp_SzCSzL);
   Matrix_Matrix_Product(Temp_SpC, Temp_SmL, Temp_SpCSmL);
   Matrix_Matrix_Product(Temp_SmC, Temp_SpL, Temp_SmCSpL);
   
   CRS Temp;
   Matrix_Matrix_Sum(Temp_SpCSmL, Temp_SmCSpL, Temp);
   Matrix_Matrix_Sum(Temp, Temp_SzCSzL, M);
   
   Matrix_Constant_Multiplication(M, coeef, 1);
   
}
