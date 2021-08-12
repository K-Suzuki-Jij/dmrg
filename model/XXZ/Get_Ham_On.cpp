#include "Model_1D_XXZ.hpp"

void Model_1D_XXZ::Get_Ham_On(CRS &M) {
   
   Check_Parameters();
   
   CRS Mag,Aniso;
   
   Model_1D_XXZ::Get_Sz_On(Mag  , h_z);
   Model_1D_XXZ::Get_SzSz_On(Aniso, D_z);
   Matrix_Matrix_Sum(Mag, Aniso, M);

}
