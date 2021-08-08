//
//  Created by Kohei Suzuki on 2020/12/23.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_SmL_On(CRS &M, int target_lspin_orbit, double coeef) {

   CRS Temp_SpL;
   
   Get_SpL_On(Temp_SpL, target_lspin_orbit, coeef);
   Matrix_Transpose(Temp_SpL, M);

}
