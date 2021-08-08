//
//  Created by Kohei Suzuki on 2020/12/23.
//
#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_SxL_On(CRS &M, int target_lspin_orbit, double coeef) {
   
   CRS Temp_SpL, Temp_SmL;
   
   Get_SpL_On(Temp_SpL, target_lspin_orbit, 0.5*coeef);
   Get_SmL_On(Temp_SmL, target_lspin_orbit, 0.5*coeef);
   Matrix_Matrix_Sum(Temp_SmL, Temp_SpL, M);
   
}

