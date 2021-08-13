//
//  main.cpp
//  1D_XXZ_ED
//
//  Created by Kohei Suzuki on 2020/04/12.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

int main(int argc, const char * argv[]) {
   
   Model_1D_XXZ Model;
   Diag_Param   Param;
   ////////////////////////////////////////////////
   Model.zero_precision = std::pow(10,-15);
   Model.p_thread = 4;
   ////////////////////////////////////////////////
   Model.BC          = "OBC";
   Model.system_size = 18;
   Model.tot_sz      = 0;
   Model.spin        = 2;
   ////////////////////////////////////////////////
   Model.J_xy = +1.0;
   Model.J_z  = +1.0;
   Model.D_z  = -0.0;
   Model.h_z  = -0.0;
   Model.site_cf_ref = 0;
   ////////////////////////////////////////////////
   Param.diag_num      = 1;
   Param.Diag_Method   = "Lanczos_Slow";
   Param.diag_acc      = pow(10,-14);
   Param.diag_min_step = 1;
   Param.diag_max_step = 800;
   ////////////////////////////////////////////////
   Param.cg_acc      = pow(10,-8);
   Param.cg_max_step = 1000;
   ////////////////////////////////////////////////
   Param.ii_Method   = "CG";
   Param.ii_acc      = pow(10,-7);
   Param.ii_diag_add = pow(10,-11);
   Param.ii_max_step = 3;
   ////////////////////////////////////////////////

   Model.dim_onsite = Model.Find_Dim_Onsite();
   Model.dim_target = Model.Find_Dim_Target(Model.tot_sz);
   
   Model.Check_Parameters();
   Model.Set_Onsite_Op();

         
   ED_Time Time;
   Time.total = omp_get_wtime();
   
   std::vector<long> Bases;
   Get_Bases(Model, Bases, Model.tot_sz, Time.bases);
   
   CRS Ham;
   Get_Ham(Model, Bases, Ham, Time.ham);
   
   std::vector<double> Vec(Model.dim_target, 0);
   Diagonalize(Ham, Vec, Param, Model, Time.diag);
   
   Onsite_Expectation_Values(Vec, Bases, Model, Time.onsite_exp);
   
   Intersite_Expectation_Values(Vec, Bases, Model, Time.intersite_exp);
   
   Time.total = omp_get_wtime() - Time.total;
   
   Model.Output_Parameters("log.txt");
   ED_Output_Time(Time, "log.txt");
   
   return 0;
}
