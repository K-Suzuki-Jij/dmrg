//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_SzC_On(CRS &M, int target_ele_orbit, double coeef) {

   CRS Temp_NC_Up, Temp_NC_Down;
   
   Get_NC_Up_On  (Temp_NC_Up  , target_ele_orbit, +0.5*coeef);
   Get_NC_Down_On(Temp_NC_Down, target_ele_orbit, -0.5*coeef);
   Matrix_Matrix_Sum(Temp_NC_Up, Temp_NC_Down, M);
   
}
