
#include "SML.hpp"
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Get_iSy_On(CRS &M, double coeef) {
   
   CRS Sp, Sm;
   
   Get_Sp_On(Sp, +0.5*coeef);
   Get_Sm_On(Sm, -0.5*coeef);

   Matrix_Matrix_Sum(Sp, Sm, M);
   
}
