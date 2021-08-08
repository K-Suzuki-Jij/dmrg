//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_NC_Up_On(CRS &M, int target_ele_orbit, double coeef) {

   CRS Temp_CUp, Temp_CUp_D;
   
   Get_CUp_On  (Temp_CUp  , target_ele_orbit, 1.0);
   Get_CUp_D_On(Temp_CUp_D, target_ele_orbit, 1.0);
   Matrix_Matrix_Product(Temp_CUp_D, Temp_CUp, M);
   Matrix_Constant_Multiplication(M, coeef, 1);

}
