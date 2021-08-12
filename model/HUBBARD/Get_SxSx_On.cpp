#include "SML.hpp"
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Get_SxSx_On(CRS &M, double coeef) {
   
    CRS Temp_Sx_On;
    Get_Sx_On(Temp_Sx_On, 1.0);
    Matrix_Matrix_Product(Temp_Sx_On, Temp_Sx_On, M);
    Matrix_Constant_Multiplication(M, coeef, 1);
    
}
