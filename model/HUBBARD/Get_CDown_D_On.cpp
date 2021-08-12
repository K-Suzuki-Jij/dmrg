#include "SML.hpp"
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Get_CDown_D_On(CRS &M, double coeef) {
   
    CRS Temp_CDown_On;
    Get_CDown_On(Temp_CDown_On, coeef);
    Matrix_Transpose(Temp_CDown_On, M);
    
}
