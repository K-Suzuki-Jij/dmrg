#include "SML.hpp"
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Get_CUp_D_On(CRS &M, double coeef) {
   
    CRS Temp_CUp_On;
    Get_CUp_On(Temp_CUp_On, coeef);
    Matrix_Transpose(Temp_CUp_On, M);
    
}
