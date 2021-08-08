//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_CUp_D_On(CRS &M, int target_ele_orbit, double coeef) {

   CRS Temp_CUp;
   
   Get_CUp_On(Temp_CUp, target_ele_orbit, coeef);
   Matrix_Transpose(Temp_CUp, M);

}
