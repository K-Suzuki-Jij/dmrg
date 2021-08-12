
#include "SML.hpp"
#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Get_SySy_On(CRS &M, double coeef) {
      
   CRS iSy;
   
   Get_iSy_On(iSy, 1.0);
   
   Matrix_Matrix_Product(iSy, iSy, M);
   
   for (int i = 0; i < M.row_dim; i++) {
      for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
         M.Val[j] = -1.0*coeef*M.Val[j];
      }
   }
   
}
