//
//  Allocate_Block.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Allocate_Block(Block_Operator &Block, Model_1D_HUBBARD &Model) {
   
   Block.Ham.resize(Model.system_size);
   Block.Sz_RE.resize(Model.system_size);
   Block.Sp_RE.resize(Model.system_size);
   Block.Sm_RE.resize(Model.system_size);
   Block.CUp_RE.resize(Model.system_size);
   Block.CDown_RE.resize(Model.system_size);
   Block.CUp_D_RE.resize(Model.system_size);
   Block.CDown_D_RE.resize(Model.system_size);
   Block.NC_RE.resize(Model.system_size);
   
   if (Model.BC == "PBC") {
      Block.Sz_LE.resize(Model.system_size);
      Block.Sp_LE.resize(Model.system_size);
      Block.Sm_LE.resize(Model.system_size);
      Block.CUp_LE.resize(Model.system_size);
      Block.CDown_LE.resize(Model.system_size);
      Block.CUp_D_LE.resize(Model.system_size);
      Block.CDown_D_LE.resize(Model.system_size);
      Block.NC_LE.resize(Model.system_size);
   }
   
}
