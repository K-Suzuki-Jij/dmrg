//
//  Allocate_Basis.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/11/18.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Allocate_Basis(DMRG_Basis_Stored &Basis, Model_1D_AKLM &Model) {
   
   Basis.LL_LLLR_Stored.resize(Model.system_size);
   Basis.LR_LLLR_Stored.resize(Model.system_size);

   Basis.Inv_LLLR_Stored.resize(Model.system_size);
   
   Basis.Trans_Mat_Stored.resize(Model.system_size);

   Basis.QN1_LL_LL_Stored.resize(Model.system_size);
   Basis.QN2_LL_LL_Stored.resize(Model.system_size);
   
   Basis.QN1_LL_LL_Stored[0].resize(Model.dim_onsite);
   Basis.QN2_LL_LL_Stored[0].resize(Model.dim_onsite);
   for (int basis = 0; basis < Model.dim_onsite; basis++) {
      Basis.QN1_LL_LL_Stored[0][basis] = Model.Find_Site_Ele(basis);
      Basis.QN2_LL_LL_Stored[0][basis] = Model.Find_Site_Sz(basis);
   }
   
}
