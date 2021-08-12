#include <cmath>
#include "Model_1D_TAKLM.hpp"

void Model_1D_TAKLM::Set_Onsite_Op() {
   
   dim_lspin  = Find_Dim_Lspin();
   dim_ele    = Find_Dim_Ele();
   dim_onsite = Find_Dim_Onsite();
   num_of_qn  = 3;
   
   Get_Ham_On       (Ham_On);
   Get_SC_1SL_On    (SC_1SL_On    , 1.0);
   Get_SC_2SL_On    (SC_2SL_On    , 1.0);
   Get_SxL_On       (SxL_On       , 1.0);
   Get_iSyL_On      (iSyL_On      , 1.0);
   Get_SzL_On       (SzL_On       , 1.0);
   Get_SpL_On       (SpL_On       , 1.0);
   Get_SmL_On       (SmL_On       , 1.0);
   Get_SxLSxL_On    (SxLSxL_On    , 1.0);
   Get_SyLSyL_On    (SyLSyL_On    , 1.0);
   Get_SzLSzL_On    (SzLSzL_On    , 1.0);
   Get_SxC_1_On     (SxC_1_On     , 1.0);
   Get_iSyC_1_On    (iSyC_1_On    , 1.0);
   Get_SzC_1_On     (SzC_1_On     , 1.0);
   Get_SpC_1_On     (SpC_1_On     , 1.0);
   Get_SmC_1_On     (SmC_1_On     , 1.0);
   Get_SxC_2_On     (SxC_2_On     , 1.0);
   Get_iSyC_2_On    (iSyC_2_On    , 1.0);
   Get_SzC_2_On     (SzC_2_On     , 1.0);
   Get_SpC_2_On     (SpC_2_On     , 1.0);
   Get_SmC_2_On     (SmC_2_On     , 1.0);
   Get_SxC_1SxC_1_On(SxC_1SxC_1_On, 1.0);
   Get_SyC_1SyC_1_On(SyC_1SyC_1_On, 1.0);
   Get_SzC_1SzC_1_On(SzC_1SzC_1_On, 1.0);
   Get_SxC_2SxC_2_On(SxC_2SxC_2_On, 1.0);
   Get_SyC_2SyC_2_On(SyC_2SyC_2_On, 1.0);
   Get_SzC_2SzC_2_On(SzC_2SzC_2_On, 1.0);
   Get_CUp_1_On     (CUp_1_On     , 1.0);
   Get_CUp_1_D_On   (CUp_1_D_On   , 1.0);
   Get_CDown_1_On   (CDown_1_On   , 1.0);
   Get_CDown_1_D_On (CDown_1_D_On , 1.0);
   Get_CUp_2_On     (CUp_2_On     , 1.0);
   Get_CUp_2_D_On   (CUp_2_D_On   , 1.0);
   Get_CDown_2_On   (CDown_2_On   , 1.0);
   Get_CDown_2_D_On (CDown_2_D_On , 1.0);
   Get_NCUp_1_On    (NCUp_1_On    , 1.0);
   Get_NCDown_1_On  (NCDown_1_On  , 1.0);
   Get_NC_1_On      (NC_1_On      , 1.0);
   Get_NCUp_2_On    (NCUp_2_On    , 1.0);
   Get_NCDown_2_On  (NCDown_2_On  , 1.0);
   Get_NC_2_On      (NC_2_On      , 1.0);
   Get_NC_On        (NC_On        , 1.0);
   
   Site_Constant.resize(system_size);
   for (int i = 0; i < system_size; i++) {
      Site_Constant[i] = std::pow(Find_Dim_Onsite(), i);
   }
   
}
