//
//  Make_Elem_Ham_LLLRRRRL.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/06/02.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Elem_Ham_LLLRRRRL(DMRG_Onsite_Basis &Basis_Onsite, DMRG_A_Basis_Set &A_Basis, DMRG_Basis &Basis, Block_Operator &System, Block_Operator &Enviro, DMRG_Block_Hamiltonian &Ham, Model_1D_XXZ &Model) {

   std::vector<int> Dummy;//This wont be used
   
   DMRG_Make_Elem_LLLR_Onsite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLR, Basis.LLLRRRRL.Inv, Ham.LLLR, 1.0);
   DMRG_Make_Elem_RRRL_Onsite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.RRRL, Basis.LLLRRRRL.Inv, Ham.RRRL, 1.0);
   DMRG_Make_Elem_LRRL_Onsite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LRRL, Basis.LLLRRRRL.Inv, Ham.LRRL, 1.0, Dummy, Dummy, "No");

   /*
   int LL_site = A_Basis.LL_site;
   int RR_site = A_Basis.RR_site;

   //Onsite Ham
   DMRG_Make_Elem_LL_Onsite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, System.Ham[LL_site], 1.0);
   DMRG_Make_Elem_LR_Onsite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Model.Ham_On, 1.0);
   DMRG_Make_Elem_RR_Onsite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Enviro.Ham[RR_site], 1.0);
   DMRG_Make_Elem_RL_Onsite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Model.Ham_On, 1.0);

   //Interaction LL-LR
   DMRG_Make_Elem_LLLR_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, System.Sp_RE[LL_site], Model.Sm_On, 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_LLLR_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, System.Sm_RE[LL_site], Model.Sp_On, 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_LLLR_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, System.Sz_RE[LL_site], Model.Sz_On,     Model.J_z , Dummy, "No");
   
   //Interaction RL-RR
   DMRG_Make_Elem_RRRL_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Model.Sp_On, Enviro.Sm_RE[RR_site], 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_RRRL_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Model.Sm_On, Enviro.Sp_RE[RR_site], 0.5*Model.J_xy, Dummy, "No");
   DMRG_Make_Elem_RRRL_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Model.Sz_On, Enviro.Sz_RE[RR_site],     Model.J_z , Dummy, "No");

   if (A_Basis.c_obc) {
      //Interaction LR-RL
      DMRG_Make_Elem_LRRL_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Model.Sp_On, Model.Sm_On, 0.5*Model.J_xy, Dummy, Dummy, "No");
      DMRG_Make_Elem_LRRL_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Model.Sm_On, Model.Sp_On, 0.5*Model.J_xy, Dummy, Dummy, "No");
      DMRG_Make_Elem_LRRL_Intersite_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv, Model.Sz_On, Model.Sz_On,     Model.J_z , Dummy, Dummy, "No");
   }
   else {
      std::cout << "Error in Make_Elem_Ham_LLLRRRRL" << std::endl;
      exit(1);
   }
   */
   
   
   DMRG_Make_Elem_Zero_LLLRRRRL(Basis_Onsite, A_Basis, Basis.LLLRRRRL.Inv);
   
   DMRG_Clear_Check_Basis(A_Basis, Basis.LLLRRRRL.Inv);
   
}
