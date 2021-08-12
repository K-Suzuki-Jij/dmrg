//
//  Find_GS_QN.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Find_GS_QN(DMRG_Quntum_Number &QN, DMRG_Block_Composition &Block_Compo, Model_1D_AKLM_TVF &Model) {
   
   QN.qn1_LLLRRRRL = (double)Model.tot_ele*(Block_Compo.LL_site + Block_Compo.RR_site + 4)/Model.system_size;
   QN.qn2_LLLRRRRL = Model.tot_parity;
   
}
