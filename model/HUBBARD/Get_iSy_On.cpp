#include "SML.hpp"
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Get_iSy_On(CRS &M, double coeef) {
    
    Check_Parameters();
    CRS Temp_SpC_On, Temp_SmC_On;
    Get_Sp_On(Temp_SpC_On, coeef*0.5);
    Get_Sm_On(Temp_SmC_On, -1.0*coeef*0.5);
    Matrix_Matrix_Sum(Temp_SpC_On, Temp_SmC_On, M);
    
}
