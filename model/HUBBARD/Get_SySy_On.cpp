#include "SML.hpp"
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Get_SySy_On(CRS &M, double coeef) {
   
    CRS Temp_iSy_On;
    Get_iSy_On(Temp_iSy_On, 1.0);
    Matrix_Matrix_Product(Temp_iSy_On, Temp_iSy_On, M);
    Matrix_Constant_Multiplication(M, -coeef, 1);
    
}
