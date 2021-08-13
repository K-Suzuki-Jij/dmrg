//
//  Get_Bases.cpp
//  1D_XXZ_ED
//
//  Created by Kohei Suzuki on 2020/04/12.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Bases(Model_1D_XXZ &Model, std::vector<long> &Bases, int target_sz, double &time) {
   
   double start = omp_get_wtime();
   
   std::vector<std::vector<int>> Basis_Seed;
   Partition_Integer((Model.system_size*Model.spin - target_sz)/2, Model.spin, Basis_Seed);
      
   std::vector<std::vector<long>> Temp_Bases(Model.p_thread);
      
   for (long i = 0; i < Basis_Seed.size(); i++) {
      int condition1 = (0 < Basis_Seed[i].size() && Basis_Seed[i].size() <= Model.system_size);
      int condition2 = (Basis_Seed[i].size() == 0 && (Model.system_size*Model.spin - target_sz)/2 == 0);
      if (condition1 || condition2) {
         
         for (long j = Basis_Seed[i].size(); j < Model.system_size; j++) {
            Basis_Seed[i].push_back(0);
         }
         
         std::vector<std::vector<int>> Temp_Basis_Seed(Model.p_thread);
         
         long basis_size = Multinomial_Coefficient(Basis_Seed[i]);
         
#pragma omp parallel num_threads (Model.p_thread)
         {
            int  thread_num = omp_get_thread_num();
            long loop_begin = thread_num*basis_size/Model.p_thread;
            long loop_end   = (thread_num + 1)*basis_size/Model.p_thread;
            
            Temp_Basis_Seed[thread_num] = Basis_Seed[i];
            Find_Nth_Permutation(Temp_Basis_Seed[thread_num], loop_begin);
            
            for (long j = loop_begin; j < loop_end; j++) {
               long basis = 0;
               for (long k = 0; k < Temp_Basis_Seed[thread_num].size(); k++) {
                  basis += Temp_Basis_Seed[thread_num][k]*Model.Site_Constant[k];
               }
               Temp_Bases[thread_num].push_back(basis);
               next_permutation(Temp_Basis_Seed[thread_num].begin(), Temp_Basis_Seed[thread_num].end());
            }
         }
      }
   }

   int dim_target = Model.Find_Dim_Target(target_sz);
   Bases.reserve(dim_target);
   
   for (int i = 0; i < Model.p_thread; i++) {
      for (long j = 0; j < Temp_Bases[i].size(); j++) {
         Bases.push_back(Temp_Bases[i][j]);
      }
   }
   
   std::sort(Bases.begin(), Bases.end());
   
   if (dim_target != Bases.size()) {
      printf("Error in Get_Bases\n");
      printf("Model.dim_target=%d, Bases.size()=%lu\n", dim_target, Bases.size());
      exit(1);
   }
   
   time = omp_get_wtime() - start;
   
}
