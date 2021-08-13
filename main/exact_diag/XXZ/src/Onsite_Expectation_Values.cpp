//
//  Expectation_Values.cpp
//  1D_XXZ_ED
//
//  Created by Kohei Suzuki on 2020/05/13.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Onsite_Expectation_Values(std::vector<double> &Eigen_Vec, std::vector<long> &Bases, Model_1D_XXZ &Model, double &time) {
   
   double start = omp_get_wtime();
   
   std::vector<double> Onsite_Val(Model.system_size, 0);
   
   ED_Expectation_Onsite(Model.Sz_On, Eigen_Vec, Bases, Onsite_Val, Model.system_size, Model.p_thread);
   Model.Output_Onsite_Values(Onsite_Val, "Sz.txt");
   Model.Output_Average_Values(Onsite_Val, "avg_Sz.txt");

   ED_Expectation_Onsite(Model.SxSx_On, Eigen_Vec, Bases, Onsite_Val, Model.system_size, Model.p_thread);
   Model.Output_Onsite_Values(Onsite_Val, "SxSx.txt");
   Model.Output_Average_Values(Onsite_Val, "avg_SxSx.txt");
   
   ED_Expectation_Onsite(Model.SySy_On, Eigen_Vec, Bases, Onsite_Val, Model.system_size, Model.p_thread);
   Model.Output_Onsite_Values(Onsite_Val, "SySy.txt");
   Model.Output_Average_Values(Onsite_Val, "avg_SySy.txt");

   ED_Expectation_Onsite(Model.SzSz_On, Eigen_Vec, Bases, Onsite_Val, Model.system_size, Model.p_thread);
   Model.Output_Onsite_Values(Onsite_Val, "SzSz.txt");
   Model.Output_Average_Values(Onsite_Val, "avg_SzSz.txt");
   
   time = omp_get_wtime() - start;

}
