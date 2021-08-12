//
//  Allocate_Block.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Allocate_Block(Block_Operator &Block, Model_1D_AKLM_TVF &Model) {
   
   Block.Ham.resize(Model.system_size);
   Block.SzL_RE.resize(Model.system_size);
   Block.SpL_RE.resize(Model.system_size);
   Block.SmL_RE.resize(Model.system_size);
   Block.SzC_RE.resize(Model.system_size);
   Block.SpC_RE.resize(Model.system_size);
   Block.SmC_RE.resize(Model.system_size);
   Block.CEven_RE.resize(Model.system_size);
   Block.COdd_RE.resize(Model.system_size);
   Block.CEven_D_RE.resize(Model.system_size);
   Block.COdd_D_RE.resize(Model.system_size);
   Block.NC_RE.resize(Model.system_size);
   
   if (Model.BC == "PBC") {
      Block.SzL_LE.resize(Model.system_size);
      Block.SpL_LE.resize(Model.system_size);
      Block.SmL_LE.resize(Model.system_size);
      Block.SzC_LE.resize(Model.system_size);
      Block.SpC_LE.resize(Model.system_size);
      Block.SmC_LE.resize(Model.system_size);
      Block.CEven_LE.resize(Model.system_size);
      Block.COdd_LE.resize(Model.system_size);
      Block.CEven_D_LE.resize(Model.system_size);
      Block.COdd_D_LE.resize(Model.system_size);
      Block.NC_LE.resize(Model.system_size);
   }
   
}
