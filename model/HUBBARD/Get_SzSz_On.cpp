#include "SML.hpp"
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Get_SzSz_On(CRS &M, double coeef) {
   
    CRS Temp_Sz_On;
    Get_Sz_On(Temp_Sz_On, 1.0);
    Matrix_Matrix_Product(Temp_Sz_On, Temp_Sz_On, M);
    Matrix_Constant_Multiplication(M, coeef, 1);
    
}
