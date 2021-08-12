#include "SML.hpp"
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Get_Ham_On(CRS &M) {
    
    CRS Sz, NCUp, NCDown, NC;
    
    Get_Sz_On(Sz, -h_z);
    Get_NCUp_On(NCUp, 1.0);
    Get_NCDown_On(NCDown, 1.0);
    Get_NC_On(NC, -mu);
    
    CRS NCUpNCDown;
    Matrix_Matrix_Product(NCUp, NCDown, NCUpNCDown);
    Matrix_Constant_Multiplication(NCUpNCDown, U, 1);
    
    CRS Temp1;
    Matrix_Matrix_Sum(Sz, NC, Temp1);
    Matrix_Matrix_Sum(Temp1, NCUpNCDown, M);
    
    Check_Symmetric_Matrix(M, zero_precision, 1);
}
