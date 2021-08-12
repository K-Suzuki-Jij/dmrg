#include <cmath>
#include "Model_1D_AKLM.hpp"

void Model_1D_AKLM::Set_Onsite_Op() {
   
   dim_lspin  = Find_Dim_Lspin();
   dim_ele    = Find_Dim_Ele();
   dim_onsite = Find_Dim_Onsite();
   num_of_qn  = 2;
      
   Get_Ham_On    (Ham_On);
   Get_SCSL_On   (SCSL_On   , 1.0);
   Get_SxL_On    (SxL_On    , 1.0);
   Get_iSyL_On   (iSyL_On   , 1.0);
   Get_SzL_On    (SzL_On    , 1.0);
   Get_SpL_On    (SpL_On    , 1.0);
   Get_SmL_On    (SmL_On    , 1.0);
   Get_SxLSxL_On (SxLSxL_On , 1.0);
   Get_SyLSyL_On (SyLSyL_On , 1.0);
   Get_SzLSzL_On (SzLSzL_On , 1.0);
   Get_SxC_On    (SxC_On    , 1.0);
   Get_iSyC_On   (iSyC_On   , 1.0);
   Get_SzC_On    (SzC_On    , 1.0);
   Get_SpC_On    (SpC_On    , 1.0);
   Get_SmC_On    (SmC_On    , 1.0);
   Get_SxCSxC_On (SxCSxC_On , 1.0);
   Get_SyCSyC_On (SyCSyC_On , 1.0);
   Get_SzCSzC_On (SzCSzC_On , 1.0);
   Get_CUp_On    (CUp_On    , 1.0);
   Get_CUp_D_On  (CUp_D_On  , 1.0);
   Get_CDown_On  (CDown_On  , 1.0);
   Get_CDown_D_On(CDown_D_On, 1.0);
   Get_NCUp_On   (NCUp_On   , 1.0);
   Get_NCDown_On (NCDown_On , 1.0);
   Get_NC_On     (NC_On     , 1.0);
   
   Site_Constant.resize(system_size);
   for (int i = 0; i < system_size; i++) {
      Site_Constant[i] = std::pow(Find_Dim_Onsite(), i);
   }
   
}
