//
//  Allocate_Block.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/05/19.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Allocate_Block(Block_Operator &Block, Model_1D_XXZ &Model) {
   
   Block.Ham.resize(Model.system_size);
   Block.Sz_RE.resize(Model.system_size);
   Block.Sp_RE.resize(Model.system_size);
   Block.Sm_RE.resize(Model.system_size);
   
   if (Model.BC == "PBC") {
      Block.Sz_LE.resize(Model.system_size);
      Block.Sp_LE.resize(Model.system_size);
      Block.Sm_LE.resize(Model.system_size);
   }
   
   Block.Ham[0]          = Model.Ham_On;
   Block.Sz_RE[0]       = Model.Sz_On;
   Block.Sp_RE[0]       = Model.Sp_On;
   Block.Sm_RE[0]       = Model.Sm_On;

   if (Model.BC == "PBC") {
      Block.Sz_LE[0]       = Model.Sz_On;
      Block.Sp_LE[0]       = Model.Sp_On;
      Block.Sm_LE[0]       = Model.Sm_On;
   }
   
}
