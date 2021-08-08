//
//  Created by Kohei Suzuki on 2020/12/16.
//

#include "SML.hpp"
#include "Model_1D_EKLM.hpp"

void Model_1D_EKLM::Get_Onsite_Operator() {
   
   num_ele_orbit   = (int)Tot_Ele.size();
   num_lspin_orbit = (int)Magnitude_2LSpin.size();
   dim_onsite      = Find_Dim_Onsite();
   
   Get_Ham_On(Ham_On);
   Get_NC_Tot_On(NC_Tot_On, 1.0);

   CUp_On    .resize(num_ele_orbit);
   CUp_D_On  .resize(num_ele_orbit);
   CDown_On  .resize(num_ele_orbit);
   CDown_D_On.resize(num_ele_orbit);
   NC_Up_On  .resize(num_ele_orbit);
   NC_Down_On.resize(num_ele_orbit);
   NC_On     .resize(num_ele_orbit);
   SpC_On    .resize(num_ele_orbit);
   SmC_On    .resize(num_ele_orbit);
   SzC_On    .resize(num_ele_orbit);
   SxC_On    .resize(num_ele_orbit);
   iSyC_On   .resize(num_ele_orbit);
   
   for (int ele_orbit = 0; ele_orbit < num_ele_orbit; ele_orbit++) {
      Get_CUp_On    (CUp_On    [ele_orbit], ele_orbit, 1.0);
      Get_CUp_D_On  (CUp_D_On  [ele_orbit], ele_orbit, 1.0);
      Get_CDown_On  (CDown_On  [ele_orbit], ele_orbit, 1.0);
      Get_CDown_D_On(CDown_D_On[ele_orbit], ele_orbit, 1.0);
      Get_NC_Up_On  (NC_Up_On  [ele_orbit], ele_orbit, 1.0);
      Get_NC_Down_On(NC_Down_On[ele_orbit], ele_orbit, 1.0);
      Get_NC_On     (NC_On     [ele_orbit], ele_orbit, 1.0);
      Get_SpC_On    (SpC_On    [ele_orbit], ele_orbit, 1.0);
      Get_SmC_On    (SmC_On    [ele_orbit], ele_orbit, 1.0);
      Get_SzC_On    (SzC_On    [ele_orbit], ele_orbit, 1.0);
      Get_SxC_On    (SxC_On    [ele_orbit], ele_orbit, 1.0);
      Get_iSyC_On   (iSyC_On   [ele_orbit], ele_orbit, 1.0);
   }
   
   SzL_On .resize(num_lspin_orbit);
   SxL_On .resize(num_lspin_orbit);
   SpL_On .resize(num_lspin_orbit);
   SmL_On .resize(num_lspin_orbit);
   iSyL_On.resize(num_lspin_orbit);

   for (int lspin_orbit = 0; lspin_orbit < num_lspin_orbit; lspin_orbit++) {
      Get_SzL_On (SzL_On [lspin_orbit], lspin_orbit, 1.0);
      Get_SxL_On (SxL_On [lspin_orbit], lspin_orbit, 1.0);
      Get_SpL_On (SpL_On [lspin_orbit], lspin_orbit, 1.0);
      Get_SmL_On (SmL_On [lspin_orbit], lspin_orbit, 1.0);
      Get_iSyL_On(iSyL_On[lspin_orbit], lspin_orbit, 1.0);
   }
   
   SzLSzL_On.resize(num_lspin_orbit);
   for (int lspin_orbit_1 = 0; lspin_orbit_1 < num_lspin_orbit; lspin_orbit_1++) {
      SzLSzL_On[lspin_orbit_1].resize(num_lspin_orbit);
      for (int lspin_orbit_2 = 0; lspin_orbit_2 < num_lspin_orbit; lspin_orbit_2++) {
         Get_SzLSzL_On(SzLSzL_On[lspin_orbit_1][lspin_orbit_2], lspin_orbit_1, lspin_orbit_2, 1.0);
      }
   }
   
   num_of_qn = num_ele_orbit + 1;
   
}
