#include "SML.hpp"

void Print_CRS(const CRS &M, std::string Name) {
   
   for (int i = 0; i < M.row_dim; i++) {
      for (long j = M.Row[i]; j < M.Row[i+1]; j++) {
         printf("%s[%-2u][%-2u]=%+.15lf\n", Name.c_str(), i, M.Col[j], M.Val[j]);
      }
   }
   
}
