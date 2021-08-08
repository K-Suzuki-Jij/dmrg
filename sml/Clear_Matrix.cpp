#include "SML.hpp"

template <class SPM> void Clear_Matrix(SPM &M) {
   
   M.Row.clear();
   M.Col.clear();
   M.Val.clear();
   M.row_dim = 0;
   M.col_dim = 0;
   
}

template void Clear_Matrix<CRS> (CRS &);
template void Clear_Matrix<CCS> (CCS &);
