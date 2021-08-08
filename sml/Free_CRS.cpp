//  Created by Kohei Suzuki on 2021/01/01.
//  Copyright © 2021 Kohei Suzuki. All rights reserved.
//

#include "SML.hpp"

void Free_CRS(CRS &M) {
   
   M.row_dim = 0;
   M.col_dim = 0;
   std::vector<double>().swap(M.Val);
   std::vector<int   >().swap(M.Col);
   std::vector<long  >().swap(M.Row);
   
}
