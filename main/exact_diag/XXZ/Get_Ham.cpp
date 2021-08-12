//
//  Get_Ham.cpp
//  1D_XXZ_ED
//
//  Created by Kohei Suzuki on 2020/04/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Ham(Model_1D_XXZ &Model, std::vector<long> &Bases, CRS &Ham, double &time) {
   
   double start = omp_get_wtime();
   
   A_Basis_Set A_Set[Model.p_thread];
   
   for (int thread_num = 0; thread_num < Model.p_thread; thread_num++) {
      A_Set[thread_num].zero_precision = Model.zero_precision;
      A_Set[thread_num].Site_Constant  = Model.Site_Constant;
      A_Set[thread_num].Local_Basis.resize(Model.system_size);
   }
   
   std::vector<long> Row_Elem_Num(Model.dim_target + 1, 0);
   
#pragma omp parallel for num_threads (Model.p_thread)
   for (int row = 0; row < Model.dim_target; row++) {
      int thread_num = omp_get_thread_num();
      
      Make_Element_Ham(Bases[row], A_Set[thread_num], Model);
      
      for (long i = 0; i < A_Set[thread_num].Basis.size(); i++) {
         long   a_basis = A_Set[thread_num].Basis[i];
         long   inv     = Binary_Search(Bases, 0, Model.dim_target, a_basis);
         if (0 <= inv && inv <= row) {
            Row_Elem_Num[row + 1]++;
         }
      }
      A_Set[thread_num].Val.clear();
      A_Set[thread_num].Basis.clear();
      A_Set[thread_num].Check_Basis.clear();
   }
   
   long tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (Model.p_thread)
   for (int row = 0; row <= Model.dim_target; row++) {
      tot_elem_num += Row_Elem_Num[row];
   }
   
   //Do not use openmp here
   for (int row = 0; row < Model.dim_target; row++) {
      Row_Elem_Num[row + 1] += Row_Elem_Num[row];
   }
   
   
   //Make Hamiltonian
   Ham.Row.resize(Model.dim_target + 1);
   Ham.Col.resize(tot_elem_num);
   Ham.Val.resize(tot_elem_num);
   
   Ham.Row[0] = 0;
#pragma omp parallel for num_threads (Model.p_thread)
   for (int row = 0; row < Model.dim_target; row++) {
      int thread_num = omp_get_thread_num();
      Make_Element_Ham(Bases[row], A_Set[thread_num], Model);
      
      for (long i = 0; i < A_Set[thread_num].Basis.size(); i++) {
         long   a_basis = A_Set[thread_num].Basis[i];
         double val     = A_Set[thread_num].Val[i];
         long   inv     = Binary_Search(Bases, 0, Model.dim_target, a_basis);
         if (0 <= inv && inv <= row) {
            Ham.Col[Row_Elem_Num[row]] = (int)inv;
            Ham.Val[Row_Elem_Num[row]] = val;
            Row_Elem_Num[row]++;
         }
      }
      Ham.Row[row + 1] = Row_Elem_Num[row];
      A_Set[thread_num].Val.clear();
      A_Set[thread_num].Basis.clear();
      A_Set[thread_num].Check_Basis.clear();
   }
   
   Ham.col_dim = Model.dim_target;
   Ham.row_dim = Model.dim_target;
   
   if (Ham.Row[Model.dim_target] != tot_elem_num) {
      printf("Error in GET_HAM\n");
      printf("Ham_Row[dim]=%lu, tot_elem_num=%lu\n",Ham.Row[Model.dim_target], tot_elem_num);
      exit(1);
   }
   if (Ham.Col.size() != tot_elem_num || Ham.Val.size() != tot_elem_num) {
      printf("Error in GET_HAM\n");
      printf("Ham_Col_size=%lu, Ham_Val_size=%lu, tot_elem_num=%lu\n", Ham.Col.size(), Ham.Val.size(), tot_elem_num);
      exit(1);
   }
   
   Sort_Col_CRS(Ham, Model.p_thread);
         
   time = omp_get_wtime() - start;
   
}
