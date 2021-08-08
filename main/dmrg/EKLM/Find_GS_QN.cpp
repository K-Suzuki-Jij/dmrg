//
//  Find_GS_QN.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/01.
//

#include "Header.hpp"

void Find_GS_QN(DMRG_Basis_LLLRRRRL &Basis, const DMRG_Block_Information &Block, const Model_1D_EKLM &Model) {
   
   int LL_site = Block.LL_site;
   int RR_site = Block.RR_site;
   Basis.QN.resize(Model.num_of_qn);
   for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
      Basis.QN[ele_orbit] = (double)Model.Tot_Ele[ele_orbit]*(LL_site + RR_site + 4)/Model.system_size;
   }
   Basis.QN[Model.num_ele_orbit] = Model.tot_sz*(LL_site + RR_site + 4)/Model.system_size;
   
}
