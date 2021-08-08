//
//  Expectation_Values.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/11.
//

#include "Header.hpp"

void Expectation_Values(const DMRG_Ground_State &GS, const DMRG_Basis_LLLRRRRL &Bases_LLLRRRRL, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, const DMRG_Block_Information &Block, const Model_1D_EKLM &Model, const DMRG_Param &Dmrg_Param) {
 
   int LL_site = Block.LL_site;
   int RR_site = Block.RR_site;
   
   int c1 = (Dmrg_Param.now_sweep == 0);
   int c2 = (Dmrg_Param.tot_sweep - Dmrg_Param.now_sweep >= 2);
   int c3 = (LL_site + RR_site + 4 != Model.system_size);
   int c4 = (LL_site != RR_site);
   
   if (c1 || c2 || c3 || c4) {
      return;
   }
   
   std::vector<CRS> Mat_System, Mat_Enviro, Mat_CF;
   std::vector<std::vector<double>> NC(Model.num_ele_orbit, std::vector<double>(Model.system_size, 0.0));
   
   for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
      std::string Output_Name = "NC" + std::to_string(ele_orbit);
      DMRG_Transform_Matrix_One(Model.NC_On[ele_orbit], Mat_System, Basis_System, Block, LL_site, Model.p_threads);
      DMRG_Transform_Matrix_One(Model.NC_On[ele_orbit], Mat_Enviro, Basis_Enviro, Block, RR_site, Model.p_threads);
      DMRG_Expectation_Onsite(NC[ele_orbit], Mat_System, Model.NC_On[ele_orbit], Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_threads);
      Output_Onsite_Values(NC[ele_orbit], Output_Name, Block, Dmrg_Param, Model);
   }
   
   
}
