//
//  Allocate_Block.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/04.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Allocate_Block(Block_Operator &Block, Model_1D_TAKLM &Model) {
   
   Block.Ham.resize(Model.system_size);
   Block.SzL_RE.resize(Model.system_size);
   Block.SpL_RE.resize(Model.system_size);
   Block.SmL_RE.resize(Model.system_size);
   Block.SzC_1_RE.resize(Model.system_size);
   Block.SpC_1_RE.resize(Model.system_size);
   Block.SmC_1_RE.resize(Model.system_size);
   Block.SzC_2_RE.resize(Model.system_size);
   Block.SpC_2_RE.resize(Model.system_size);
   Block.SmC_2_RE.resize(Model.system_size);
   Block.CUp_1_RE.resize(Model.system_size);
   Block.CDown_1_RE.resize(Model.system_size);
   Block.CUp_1_D_RE.resize(Model.system_size);
   Block.CDown_1_D_RE.resize(Model.system_size);
   Block.CUp_2_RE.resize(Model.system_size);
   Block.CDown_2_RE.resize(Model.system_size);
   Block.CUp_2_D_RE.resize(Model.system_size);
   Block.CDown_2_D_RE.resize(Model.system_size);
   Block.NC_1_RE.resize(Model.system_size);
   Block.NC_2_RE.resize(Model.system_size);
   Block.NC_RE.resize(Model.system_size);
   
   Block.Ham[0]          = Model.Ham_On;
   Block.SzL_RE[0]       = Model.SzL_On;
   Block.SpL_RE[0]       = Model.SpL_On;
   Block.SmL_RE[0]       = Model.SmL_On;
   Block.SzC_1_RE[0]     = Model.SzC_1_On;
   Block.SpC_1_RE[0]     = Model.SpC_1_On;
   Block.SmC_1_RE[0]     = Model.SmC_1_On;
   Block.SzC_2_RE[0]     = Model.SzC_2_On;
   Block.SpC_2_RE[0]     = Model.SpC_2_On;
   Block.SmC_2_RE[0]     = Model.SmC_2_On;
   Block.CUp_1_RE[0]     = Model.CUp_1_On;
   Block.CDown_1_RE[0]   = Model.CDown_1_On;
   Block.CUp_1_D_RE[0]   = Model.CUp_1_D_On;
   Block.CDown_1_D_RE[0] = Model.CDown_1_D_On;
   Block.CUp_2_RE[0]     = Model.CUp_2_On;
   Block.CDown_2_RE[0]   = Model.CDown_2_On;
   Block.CUp_2_D_RE[0]   = Model.CUp_2_D_On;
   Block.CDown_2_D_RE[0] = Model.CDown_2_D_On;
   Block.NC_1_RE[0]      = Model.NC_1_On;
   Block.NC_2_RE[0]      = Model.NC_2_On;
   Block.NC_RE[0]        = Model.NC_On;

   
   if (Model.BC == "PBC") {
      Block.SzL_LE.resize(Model.system_size);
      Block.SpL_LE.resize(Model.system_size);
      Block.SmL_LE.resize(Model.system_size);
      Block.SzC_1_LE.resize(Model.system_size);
      Block.SpC_1_LE.resize(Model.system_size);
      Block.SmC_1_LE.resize(Model.system_size);
      Block.SzC_2_LE.resize(Model.system_size);
      Block.SpC_2_LE.resize(Model.system_size);
      Block.SmC_2_LE.resize(Model.system_size);
      Block.CUp_1_LE.resize(Model.system_size);
      Block.CDown_1_LE.resize(Model.system_size);
      Block.CUp_1_D_LE.resize(Model.system_size);
      Block.CDown_1_D_LE.resize(Model.system_size);
      Block.CUp_2_LE.resize(Model.system_size);
      Block.CDown_2_LE.resize(Model.system_size);
      Block.CUp_2_D_LE.resize(Model.system_size);
      Block.CDown_2_D_LE.resize(Model.system_size);
      Block.NC_1_LE.resize(Model.system_size);
      Block.NC_2_LE.resize(Model.system_size);
      Block.NC_LE.resize(Model.system_size);
      
      Block.SzL_LE[0]       = Model.SzL_On;
      Block.SpL_LE[0]       = Model.SpL_On;
      Block.SmL_LE[0]       = Model.SmL_On;
      Block.SzC_1_LE[0]     = Model.SzC_1_On;
      Block.SpC_1_LE[0]     = Model.SpC_1_On;
      Block.SmC_1_LE[0]     = Model.SmC_1_On;
      Block.SzC_2_LE[0]     = Model.SzC_2_On;
      Block.SpC_2_LE[0]     = Model.SpC_2_On;
      Block.SmC_2_LE[0]     = Model.SmC_2_On;
      Block.CUp_1_LE[0]     = Model.CUp_1_On;
      Block.CDown_1_LE[0]   = Model.CDown_1_On;
      Block.CUp_1_D_LE[0]   = Model.CUp_1_D_On;
      Block.CDown_1_D_LE[0] = Model.CDown_1_D_On;
      Block.CUp_2_LE[0]     = Model.CUp_2_On;
      Block.CDown_2_LE[0]   = Model.CDown_2_On;
      Block.CUp_2_D_LE[0]   = Model.CUp_2_D_On;
      Block.CDown_2_D_LE[0] = Model.CDown_2_D_On;
      Block.NC_1_LE[0]      = Model.NC_1_On;
      Block.NC_2_LE[0]      = Model.NC_2_On;
      Block.NC_LE[0]        = Model.NC_On;
   }
   
   

}
