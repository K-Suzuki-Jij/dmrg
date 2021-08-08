//
//  DMRG_Get_Inv_LLLRRRRL.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/12/31.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include <vector>
#include "DMRG.hpp"

void DMRG_Get_Inv_LLLRRRRL(DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Block_Information &Block, int p_threads) {
   
   int dim_LL     = Block.dim_LL;
   int dim_RR     = Block.dim_RR;
   int dim_onsite = Block.dim_onsite;

   long whole_dim = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   Basis_LLLRRRRL.Inv.resize(whole_dim);
   
#pragma omp parallel for num_threads (p_threads)
   for (long i = 0; i < whole_dim; i++) {
      Basis_LLLRRRRL.Inv[i] = -1;
   }
   
   int base_LRRRRL  = dim_onsite*dim_RR*dim_onsite;
   int base_RRRL    = dim_RR*dim_onsite;
   int base_RL      = dim_onsite;
   int dim_LLLRRRRL = (int)Basis_LLLRRRRL.LL.size();
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      int  LL       = Basis_LLLRRRRL.LL[i];
      int  LR       = Basis_LLLRRRRL.LR[i];
      int  RR       = Basis_LLLRRRRL.RR[i];
      int  RL       = Basis_LLLRRRRL.RL[i];
      long LLLRRRRL = (long)LL*base_LRRRRL + LR*base_RRRL + RR*dim_onsite + RL;
      Basis_LLLRRRRL.Inv[LLLRRRRL] = i;
   }
   
   Basis_LLLRRRRL.base_LRRRRL = base_LRRRRL;
   Basis_LLLRRRRL.base_RRRL   = base_RRRL;
   Basis_LLLRRRRL.base_RL     = base_RL;
   
}
