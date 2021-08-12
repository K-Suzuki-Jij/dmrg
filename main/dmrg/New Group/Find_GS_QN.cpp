//
//  Find_GS_QN.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/04.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Find_GS_QN(int &qn1, int &qn2, int &qn3, DMRG_Block_Information &Block, Model_1D_TAKLM &Model) {
   
   int LL_site = Block.LL_site;
   int RR_site = Block.RR_site;
   
   qn1 = (double)Model.tot_ele_1*(LL_site + RR_site + 4)/Model.system_size;
   qn2 = (double)Model.tot_ele_2*(LL_site + RR_site + 4)/Model.system_size;
   qn3 = (double)Model.tot_sz   *(LL_site + RR_site + 4)/Model.system_size;
   
  if (qn1 > qn2) {
      if ((qn1 + qn2)%2 == 1 && qn3%2 == 0) {
         qn1 = qn1 + 1;
      }
      if ((qn1 + qn2)%2 == 0 && qn3%2 == 1) {
         qn1 = qn1 + 1;
      }
   }
   else {
      if ((qn1 + qn2)%2 == 1 && qn3%2 == 0) {
         qn2 = qn2 + 1;
      }
      if ((qn1 + qn2)%2 == 0 && qn3%2 == 1) {
         qn2 = qn2 + 1;
      }
   }
      
}
