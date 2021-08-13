//
//  Diagonalize.cpp
//  1D_XXZ_ED
//
//  Created by Kohei Suzuki on 2020/05/16.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Diagonalize(CRS &Ham, std::vector<double> &GS_Vec, Diag_Param &Param, Model_1D_XXZ &Model, double &time) {
   
   double start = omp_get_wtime();
   
   double GS_val = 0;
   Param.Mat_Type = "Sym";
   Param.Lanczos_Initial_Guess = "No";
   Param.Calc_Vec = "Yes";
   Lanczos(Ham, GS_Vec, GS_val, Param, Model.p_thread);
   Inverse_Iteration(Ham, GS_Vec, GS_val, Param, Model.p_thread);
   Model.Output_Average_Values(GS_val, "energy.txt");
   Lanczos_Ex1(Ham, GS_Vec, GS_val, Param, Model.p_thread);
   Model.Output_Average_Values(GS_val, "energy_ex1.txt");
   
   time = omp_get_wtime() - start;
   
}
