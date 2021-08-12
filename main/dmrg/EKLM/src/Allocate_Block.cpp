//
//  Allocate_Block.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/01.
//

#include "Header.hpp"

void Allocate_Block(Block_Operator &Block, const Model_1D_EKLM &Model) {
   
   int system_size          = Model.system_size;
   int num_ele_orbit        = Model.num_ele_orbit;
   int num_lspin_orbit      = Model.num_lspin_orbit;
   
   Block.Ham.resize(system_size);
   Block.Ham[0] = Model.Ham_On;
   
   Block.CUp    .resize(system_size);
   Block.CDown  .resize(system_size);
   Block.CUp_D  .resize(system_size);
   Block.CDown_D.resize(system_size);
   Block.NC_Tot .resize(system_size);
   
   for (int i = 0; i < system_size; i++) {
      Block.CUp    [i].resize(Model.Ele_Ele_t_site.size());
      Block.CDown  [i].resize(Model.Ele_Ele_t_site.size());
      Block.CUp_D  [i].resize(Model.Ele_Ele_t_site.size());
      Block.CDown_D[i].resize(Model.Ele_Ele_t_site.size());
      
      for (auto j = 0; j < Model.Ele_Ele_t_site.size(); j++) {
         Block.CUp    [i][j].resize(num_ele_orbit);
         Block.CDown  [i][j].resize(num_ele_orbit);
         Block.CUp_D  [i][j].resize(num_ele_orbit);
         Block.CDown_D[i][j].resize(num_ele_orbit);
      }
      
      Block.NC_Tot [i].resize(Model.Ele_Ele_V_site.size());

   }
   
   Block.NC_Tot[0][0] = Model.NC_Tot_On;
   
   for (int ele_orbit = 0; ele_orbit < num_ele_orbit; ele_orbit++) {
      Block.CUp    [0][0][ele_orbit] = Model.CUp_On    [ele_orbit];
      Block.CDown  [0][0][ele_orbit] = Model.CDown_On  [ele_orbit];
      Block.CUp_D  [0][0][ele_orbit] = Model.CUp_D_On  [ele_orbit];
      Block.CDown_D[0][0][ele_orbit] = Model.CDown_D_On[ele_orbit];
   }
 
   Block.SpL   .resize(system_size);
   Block.SmL   .resize(system_size);
   Block.SzL   .resize(system_size);
   
   for (int i = 0; i < system_size; i++) {
      Block.SpL[i].resize(Model.LSpin_LSpin_Jxy_site.size());
      Block.SmL[i].resize(Model.LSpin_LSpin_Jxy_site.size());
      
      for (auto j = 0; j < Model.LSpin_LSpin_Jxy_site.size(); j++) {
         Block.SpL[i][j].resize(num_lspin_orbit);
         Block.SmL[i][j].resize(num_lspin_orbit);
      }
      
      Block.SzL[i].resize(Model.LSpin_LSpin_Jz_site.size());
      
      for (auto j = 0; j < Model.LSpin_LSpin_Jz_site.size(); j++) {
         Block.SzL[i][j].resize(num_lspin_orbit);
      }
   }
   
   for (int lspin_orbit = 0; lspin_orbit < num_lspin_orbit; lspin_orbit++) {
      Block.SpL[0][0][lspin_orbit] = Model.SpL_On[lspin_orbit];
      Block.SmL[0][0][lspin_orbit] = Model.SmL_On[lspin_orbit];
      Block.SzL[0][0][lspin_orbit] = Model.SzL_On[lspin_orbit];
   }
   
}
