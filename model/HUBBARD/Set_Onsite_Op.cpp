#include <cmath>
#include "Model_1D_HUBBARD.hpp"

void Model_1D_HUBBARD::Set_Onsite_Op() {
    
    dim_onsite = Find_Dim_Onsite();
    num_of_qn  = 2;
    
    Get_Ham_On    (Ham_On    );
    Get_Sx_On     (Sx_On     , 1.0);
    Get_iSy_On    (iSy_On    , 1.0);
    Get_Sz_On     (Sz_On     , 1.0);
    Get_Sp_On     (Sp_On     , 1.0);
    Get_Sm_On     (Sm_On     , 1.0);
    Get_SxSx_On   (SxSx_On   , 1.0);
    Get_SySy_On   (SySy_On   , 1.0);
    Get_SzSz_On   (SzSz_On   , 1.0);
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
