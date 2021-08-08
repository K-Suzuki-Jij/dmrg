//
//  Created by Kohei Suzuki on 2020/12/23.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_SzLSzL_On(CRS &M, int target_lspin_orbit_1, int target_lspin_orbit_2, double coeef) {
   
   CRS Temp_SzL1, Temp_SzL2;
   Get_SzL_On(Temp_SzL1, target_lspin_orbit_1, 1.0);
   Get_SzL_On(Temp_SzL2, target_lspin_orbit_2, 1.0);
   Matrix_Matrix_Product(Temp_SzL1, Temp_SzL2, M);
   Matrix_Constant_Multiplication(M, coeef, 1);

}
