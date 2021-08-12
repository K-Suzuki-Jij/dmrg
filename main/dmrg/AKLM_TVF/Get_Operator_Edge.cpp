//
//  Get_Operator_Edge.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Operator_Edge(Block_Operator &Block, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_Stored, Model_1D_AKLM_TVF &Model) {
   
   Block.Ham[0]        = Model.Ham_On;
   Block.SzL_RE[0]     = Model.SzL_On;
   Block.SpL_RE[0]     = Model.SpL_On;
   Block.SmL_RE[0]     = Model.SmL_On;
   Block.SzC_RE[0]     = Model.SzC_On;
   Block.SpC_RE[0]     = Model.SpC_On;
   Block.SmC_RE[0]     = Model.SmC_On;
   Block.CEven_RE[0]   = Model.CEven_On;
   Block.COdd_RE[0]    = Model.COdd_On;
   Block.CEven_D_RE[0] = Model.CEven_D_On;
   Block.COdd_D_RE[0]  = Model.COdd_D_On;
   Block.NC_RE[0]      = Model.NC_On;

   if (Model.BC == "PBC") {
      Block.SzL_LE[0]     = Model.SzL_On;
      Block.SpL_LE[0]     = Model.SpL_On;
      Block.SmL_LE[0]     = Model.SmL_On;
      Block.SzC_LE[0]     = Model.SzC_On;
      Block.SpC_LE[0]     = Model.SpC_On;
      Block.SmC_LE[0]     = Model.SmC_On;
      Block.CEven_LE[0]   = Model.CEven_On;
      Block.COdd_LE[0]    = Model.COdd_On;
      Block.CEven_D_LE[0] = Model.CEven_D_On;
      Block.COdd_D_LE[0]  = Model.COdd_D_On;
      Block.NC_LE[0]      = Model.NC_On;
   }
   
   Basis_Stored.QN1_LL_LL_Stored[0].resize(Model.dim_onsite);
   Basis_Stored.QN2_LL_LL_Stored[0].resize(Model.dim_onsite);
   for (int basis = 0; basis < Model.dim_onsite; basis++) {
      Basis_Stored.QN1_LL_LL_Stored[0][basis] = Model.Find_Site_Ele(basis);
      Basis_Stored.QN2_LL_LL_Stored[0][basis] = Model.Find_Site_Parity(basis);
   }
   
   Block_Compo.dim_onsite = Model.dim_onsite;
   
}
