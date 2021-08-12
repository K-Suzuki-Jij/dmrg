#include <cmath>
#include "SML.hpp"
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Get_Sm_On(CRS &M, double coeef) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
    Check_Parameters();
    
    Free_CRS(M);
    M.row_dim = dim_onsite;
    M.col_dim = dim_onsite;

    M.Row.push_back(0);
    for (int row = 0; row < dim_onsite; row++) {
        for (int col = 0; col < dim_onsite; col++) {
            double val = 0.0;
            if (col == 1 && row == 2) {
                val = 1.0;
            }
            val *= coeef;
            if (std::fabs(val) > zero_precision) {
                M.Col.push_back(col);
                M.Val.push_back(val);
            }
        }
        M.Row.push_back(M.Col.size());
    }
}
