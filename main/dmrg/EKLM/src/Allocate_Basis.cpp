//
//  Allocate_Basis.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/01.
//

#include "Header.hpp"

void Allocate_Basis(DMRG_Basis_Stored &Basis, Model_1D_EKLM &Model) {
   
   Basis.LL_LLLR_Stored  .resize(Model.system_size);
   Basis.LR_LLLR_Stored  .resize(Model.system_size);
   Basis.Inv_LLLR_Stored .resize(Model.system_size);
   Basis.Trans_Mat_Stored.resize(Model.system_size);
   
   Basis.QN_LL_LL_Stored.resize(Model.system_size);
   Basis.QN_LL_LL_Stored[0].resize(Model.num_of_qn);
   
   for (int i = 0; i < Model.num_of_qn; i++) {
      Basis.QN_LL_LL_Stored[0][i].resize(Model.dim_onsite);
   }
   
   for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
      for (int basis = 0; basis < Model.dim_onsite; basis++) {
         Basis.QN_LL_LL_Stored[0][ele_orbit][basis] = Model.Find_Site_Ele(basis, ele_orbit);
      }
   }
   
   Basis.qn_LL_LL_stored_ele_start = 0;
   Basis.qn_LL_LL_stored_ele_end   = Model.num_ele_orbit;
   
   Basis.qn_LL_LL_stored_parity_start = 0;
   Basis.qn_LL_LL_stored_parity_end   = 0;

   for (int basis = 0; basis < Model.dim_onsite; basis++) {
      Basis.QN_LL_LL_Stored[0][Model.num_ele_orbit][basis] = Model.Find_Site_Sz(basis);
   }   
   
}
