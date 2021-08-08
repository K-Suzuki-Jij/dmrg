//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_CDown_D_On(CRS &M, int target_ele_orbit, double coeef) {

   CRS Temp_CDown;
   
   Get_CDown_On(Temp_CDown, target_ele_orbit, coeef);
   Matrix_Transpose(Temp_CDown, M);

}
