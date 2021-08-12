//
//  Get_Operator_Edge.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Get_Operator_Edge(Block_Operator &Block, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_Stored, Model_1D_HUBBARD &Model) {
   
   Block.Ham[0]        = Model.Ham_On;
   Block.Sz_RE[0]      = Model.Sz_On;
   Block.Sp_RE[0]      = Model.Sp_On;
   Block.Sm_RE[0]      = Model.Sm_On;
   Block.CUp_RE[0]     = Model.CUp_On;
   Block.CDown_RE[0]   = Model.CDown_On;
   Block.CUp_D_RE[0]   = Model.CUp_D_On;
   Block.CDown_D_RE[0] = Model.CDown_D_On;
   Block.NC_RE[0]      = Model.NC_On;

   if (Model.BC == "PBC") {
      Block.Sz_LE[0]      = Model.Sz_On;
      Block.Sp_LE[0]      = Model.Sp_On;
      Block.Sm_LE[0]      = Model.Sm_On;
      Block.CUp_LE[0]     = Model.CUp_On;
      Block.CDown_LE[0]   = Model.CDown_On;
      Block.CUp_D_LE[0]   = Model.CUp_D_On;
      Block.CDown_D_LE[0] = Model.CDown_D_On;
      Block.NC_LE[0]      = Model.NC_On;
   }
   
   Basis_Stored.QN1_LL_LL_Stored[0].resize(Model.dim_onsite);
   Basis_Stored.QN2_LL_LL_Stored[0].resize(Model.dim_onsite);
   for (int basis = 0; basis < Model.dim_onsite; basis++) {
      Basis_Stored.QN1_LL_LL_Stored[0][basis] = Model.Find_Site_Ele(basis);
      Basis_Stored.QN2_LL_LL_Stored[0][basis] = Model.Find_Site_Sz(basis);
   }
   
   Block_Compo.dim_onsite = Model.dim_onsite;
   
}
