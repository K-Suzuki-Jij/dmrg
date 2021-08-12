#include <climits>
#include "Model_1D_XXZ.hpp"

int Model_1D_XXZ::Find_Dim_Target(int target_sz) {
   
   int  tot_spin = system_size*spin;
   long Dim[system_size][tot_spin + 1];
   
   if (std::abs(target_sz) > tot_spin) {
      return 0;
   }
   
   //Initialize
   for (int site = 0; site < system_size; site++) {
      for (int s = 0; s <= tot_spin; s++) {
         Dim[site][s] = 0;
      }
   }
   
   for (int s = -spin; s <= (int)spin; s+=2) {
      Dim[0][(s + spin)/2] = 1;
   }
   
   for (int site = 1; site < system_size; site++) {
      for (int s = -spin; s <= (int)spin; s+=2) {
         for (int s_prev = -spin*site; s_prev <= (int)(spin*site); s_prev+=2) {
            Dim[site][(s + s_prev + spin*(site + 1))/2] += Dim[site - 1][(s_prev + spin*site)/2];
            if (Dim[site][(s + s_prev + spin*(site + 1))/2] > LONG_MAX) {
               printf("Too large Hilbert space\n");
               exit(1);
            }
         }
      }
   }
   
   if (Dim[system_size - 1][(target_sz + tot_spin)/2] > INT_MAX) {
      printf("Too large Hilbert space\n");
      exit(1);
   }
   
   return (int)Dim[system_size - 1][(target_sz + tot_spin)/2];
}
