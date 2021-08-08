//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_SxC_On(CRS &M, int target_ele_orbit, double coeef) {

   CRS Temp_SpC, Temp_SmC;
   
   Get_SpC_On(Temp_SpC, target_ele_orbit, 0.5*coeef);
   Get_SmC_On(Temp_SmC, target_ele_orbit, 0.5*coeef);
   Matrix_Matrix_Sum(Temp_SmC, Temp_SpC, M);
   
}
