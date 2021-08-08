#include "SML.hpp"

void Print_CCS(const CCS &M, std::string Name) {
   
   for (int i = 0; i < M.row_dim; i++) {
      for (int j = 0; j < M.col_dim; j++) {
         for (long k = M.Col[j]; k < M.Col[j+1]; k++) {
            if (M.Row[k] == i) {
               printf("%s[%-2u][%-2u]=%+.15lf\n", Name.c_str(), M.Row[k], j, M.Val[k]);
            }
         }
      }
   }
   
}
