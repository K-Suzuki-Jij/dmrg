//
//  Intersite_Expectation_Values.cpp
//  1D_XXZ_ED
//
//  Created by Kohei Suzuki on 2020/05/17.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Intersite_Expectation_Values(std::vector<double> &Eigen_Vec, std::vector<long> &Bases, Model_1D_XXZ &Model, double &time) {
   
   double start = omp_get_wtime();
   
   std::vector<double> Intersite_Val(Model.system_size - Model.site_cf_ref, 0);
   
   ED_Expectation_Intersite_Q1(Model.Sz_On, Eigen_Vec, Bases, Bases, Intersite_Val, Model.site_cf_ref, Model.system_size, Model.p_thread);
   Model.Output_Intersite_Values(Intersite_Val, "Sz_CF.txt");
   
   std::vector<long> Bases_Sz_P, Bases_Sz_M;

   Get_Bases(Model, Bases_Sz_M, Model.tot_sz - 2, time);
   Get_Bases(Model, Bases_Sz_P, Model.tot_sz + 2, time);
   
   ED_Expectation_Intersite_Q1(Model.Sx_On, Eigen_Vec, Bases, Bases_Sz_P, Bases_Sz_M, Intersite_Val, Model.site_cf_ref, Model.system_size, Model.p_thread);
   Model.Output_Intersite_Values(Intersite_Val, "Sx_CF.txt");
   
   ED_Expectation_Intersite_Q1(Model.iSy_On, Eigen_Vec, Bases, Bases_Sz_P, Bases_Sz_M, Intersite_Val, Model.site_cf_ref, Model.system_size, Model.p_thread);
   Model.Output_Intersite_Values(Intersite_Val, "Sy_CF.txt");
   
   time = omp_get_wtime() - start;

}

