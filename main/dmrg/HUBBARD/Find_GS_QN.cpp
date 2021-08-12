//
//  Find_GS_QN.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Find_GS_QN(DMRG_Quntum_Number &QN, DMRG_Block_Composition &Block_Compo, Model_1D_HUBBARD &Model) {
   
   QN.qn1_LLLRRRRL = (double)Model.tot_ele*(Block_Compo.LL_site + Block_Compo.RR_site + 4)/Model.system_size;
   QN.qn2_LLLRRRRL = (double)Model.tot_sz *(Block_Compo.LL_site + Block_Compo.RR_site + 4)/Model.system_size;
   
   if (QN.qn1_LLLRRRRL%2 == 1 && QN.qn2_LLLRRRRL%2 == 0) {
      QN.qn1_LLLRRRRL = QN.qn1_LLLRRRRL + 1;
   }
   if (QN.qn1_LLLRRRRL%2 == 0 && QN.qn2_LLLRRRRL%2 == 1) {
      QN.qn1_LLLRRRRL = QN.qn1_LLLRRRRL + 1;
   }
   
   if (QN.qn1_LLLRRRRL == 0) {
      QN.qn1_LLLRRRRL = QN.qn1_LLLRRRRL + 1;
      QN.qn2_LLLRRRRL = QN.qn2_LLLRRRRL + 1;
   }
   
}
