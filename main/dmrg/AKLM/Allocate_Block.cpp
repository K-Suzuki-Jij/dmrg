//
//  Allocate_Block.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/07/08.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Allocate_Block(Block_Operator &Block, Model_1D_AKLM &Model) {
   
   Block.Ham.resize(Model.system_size);
   Block.SzL_RE.resize(Model.system_size);
   Block.SpL_RE.resize(Model.system_size);
   Block.SmL_RE.resize(Model.system_size);
   Block.SzC_RE.resize(Model.system_size);
   Block.SpC_RE.resize(Model.system_size);
   Block.SmC_RE.resize(Model.system_size);
   Block.CUp_RE.resize(Model.system_size);
   Block.CDown_RE.resize(Model.system_size);
   Block.CUp_D_RE.resize(Model.system_size);
   Block.CDown_D_RE.resize(Model.system_size);
   Block.NC_RE.resize(Model.system_size);
   
   Block.Ham[0]        = Model.Ham_On;
   Block.SzL_RE[0]     = Model.SzL_On;
   Block.SpL_RE[0]     = Model.SpL_On;
   Block.SmL_RE[0]     = Model.SmL_On;
   Block.SzC_RE[0]     = Model.SzC_On;
   Block.SpC_RE[0]     = Model.SpC_On;
   Block.SmC_RE[0]     = Model.SmC_On;
   Block.CUp_RE[0]     = Model.CUp_On;
   Block.CDown_RE[0]   = Model.CDown_On;
   Block.CUp_D_RE[0]   = Model.CUp_D_On;
   Block.CDown_D_RE[0] = Model.CDown_D_On;
   Block.NC_RE[0]      = Model.NC_On;
   
   if (Model.BC == "PBC") {
      Block.SzL_LE.resize(Model.system_size);
      Block.SpL_LE.resize(Model.system_size);
      Block.SmL_LE.resize(Model.system_size);
      Block.SzC_LE.resize(Model.system_size);
      Block.SpC_LE.resize(Model.system_size);
      Block.SmC_LE.resize(Model.system_size);
      Block.CUp_LE.resize(Model.system_size);
      Block.CDown_LE.resize(Model.system_size);
      Block.CUp_D_LE.resize(Model.system_size);
      Block.CDown_D_LE.resize(Model.system_size);
      Block.NC_LE.resize(Model.system_size);
      
      Block.SzL_LE[0]     = Model.SzL_On;
      Block.SpL_LE[0]     = Model.SpL_On;
      Block.SmL_LE[0]     = Model.SmL_On;
      Block.SzC_LE[0]     = Model.SzC_On;
      Block.SpC_LE[0]     = Model.SpC_On;
      Block.SmC_LE[0]     = Model.SmC_On;
      Block.CUp_LE[0]     = Model.CUp_On;
      Block.CDown_LE[0]   = Model.CDown_On;
      Block.CUp_D_LE[0]   = Model.CUp_D_On;
      Block.CDown_D_LE[0] = Model.CDown_D_On;
      Block.NC_LE[0]      = Model.NC_On;
   }
   
}
