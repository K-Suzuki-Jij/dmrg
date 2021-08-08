#include "SML.hpp"

void Matrix_Transpose(const CRS &M, CRS &M_Out) {
   
   int row_dim = M.row_dim;
   int col_dim = M.col_dim;
   std::vector<int> Row_Count(row_dim, 0);
   
   Clear_Matrix(M_Out);
   
   M_Out.Row.push_back(0);
   
   for (int i = 0; i < col_dim; i++) {
      for (int j = 0; j < row_dim; j++) {
         if ((M.Row[j] + Row_Count[j] < M.Row[j+1]) && (M.Col[M.Row[j] + Row_Count[j]] == i)) {
            M_Out.Val.push_back(M.Val[M.Row[j] + Row_Count[j]]);
            M_Out.Col.push_back(j);
            Row_Count[j]++;
         }
      }
      M_Out.Row.push_back(M_Out.Col.size());
   }
   
   M_Out.row_dim = col_dim;
   M_Out.col_dim = row_dim;
   
}

void Matrix_Transpose(const CCS &M, CCS &M_Out) {
   
   int row_dim = M.row_dim;
   int col_dim = M.col_dim;
   
   std::vector<int> Col_Count(col_dim, 0);
   
   Clear_Matrix(M_Out);
   
   M_Out.Col.push_back(0);
   
   for (int i = 0; i < row_dim; i++) {
      for (int j = 0; j < col_dim; j++) {
         if (M.Col[j] + Col_Count[j] < M.Col[j+1] && M.Row[M.Col[j] + Col_Count[j]] == i) {
            M_Out.Val.push_back(M.Val[M.Col[j] + Col_Count[j]]);
            M_Out.Row.push_back(j);
            Col_Count[j]++;
         }
      }
      M_Out.Col.push_back(M_Out.Row.size());
   }
   
   M_Out.row_dim = col_dim;
   M_Out.col_dim = row_dim;
   
}
