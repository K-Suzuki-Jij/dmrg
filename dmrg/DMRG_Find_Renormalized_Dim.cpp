#include <iostream>
#include <cmath>
#include "DMRG.hpp"

int DMRG_Find_Renormalized_Dim(std::vector<std::vector<double>> &Value_DM, std::vector<int> &Block_Sorted, std::vector<int> &Basis_Sorted, int max_dim_system) {
   
   int dim_LLLR   = (int)Block_Sorted.size();
   int dim_renorm = 0;
   double zero_precision = std::pow(10,-15);
   
   if (dim_LLLR <= max_dim_system) {
      for (int i = dim_LLLR; i >= 1; i--) {
         double val = Value_DM[Block_Sorted[i - 1]][Basis_Sorted[i - 1]];
         if (val > zero_precision) {
            dim_renorm = i;
            break;
         }
      }
   }
   else {
      for (int i = max_dim_system; i >= 1; i--) {
         double val1 = Value_DM[Block_Sorted[i - 1]][Basis_Sorted[i - 1]];
         double val2 = Value_DM[Block_Sorted[i    ]][Basis_Sorted[i    ]];
         if (std::abs(val1 - val2) > zero_precision && val1 > zero_precision) {
            dim_renorm = i;
            break;
         }
      }
   }
   
   if (dim_renorm <= 0) {
      std::cout << "Error in DMRG_Find_Renormalized_Dim" << std::endl;
      exit(1);
   }
   
   return dim_renorm;
   
}

