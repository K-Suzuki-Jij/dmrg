//
//  Make_Elem_Ham.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/04.
//

#include "Header.hpp"

void Make_Elem_Ham(std::string Mat_Type, const DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv, const Block_Operator &System, const Block_Operator &Enviro, const std::vector<int> &Ele, const Model_1D_EKLM &Model) {
   
   int LL_site = A_Basis.LL_site;
   int RR_site = A_Basis.RR_site;
   
   if (Mat_Type == "LLLR") {
      
      //Onsite Ham
      DMRG_Make_Elem_Onsite_LLLR_RRRL("LL", Basis.LL, Basis.LR, A_Basis, Inv, System.Ham[LL_site], 1.0);
      DMRG_Make_Elem_Onsite_LLLR_RRRL("LR", Basis.LR, Basis.LL, A_Basis, Inv, Model.Ham_On       , 1.0);
      
      //Interaction: Hopping
      for (int i = 0; i < Model.Ele_Ele_t_site.size() && LL_site - i >= 0; i++) {
         for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
            DMRG_Make_Elem_Intersite_Block("LL_LR", Basis.LL, Basis.LR, A_Basis, Inv, System.CUp_D  [LL_site][i][ele_orbit], Model.CUp_On    [ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LL_LR", Basis.LL, Basis.LR, A_Basis, Inv, System.CUp    [LL_site][i][ele_orbit], Model.CUp_D_On  [ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LL_LR", Basis.LL, Basis.LR, A_Basis, Inv, System.CDown_D[LL_site][i][ele_orbit], Model.CDown_On  [ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LL_LR", Basis.LL, Basis.LR, A_Basis, Inv, System.CDown  [LL_site][i][ele_orbit], Model.CDown_D_On[ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
         }
      }
      
      //Interaction: Intersite Coulomb
      for (int i = 0; i < Model.Ele_Ele_V_site.size() && LL_site - i >= 0; i++) {
         DMRG_Make_Elem_Intersite_Block("LL_LR", Basis.LL, Basis.LR, A_Basis, Inv, System.NC_Tot[LL_site][i], Model.NC_Tot_On, Model.Ele_Ele_V_site[i], Ele, "No");
      }
      
      //Interaction: Local Spin
      for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
         for (int i = 0; i < Model.LSpin_LSpin_Jxy_site.size() && LL_site - i >= 0; i++) {
            DMRG_Make_Elem_Intersite_Block("LL_LR", Basis.LL, Basis.LR, A_Basis, Inv, System.SpL[LL_site][i][lspin_orbit], Model.SmL_On[lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
            DMRG_Make_Elem_Intersite_Block("LL_LR", Basis.LL, Basis.LR, A_Basis, Inv, System.SmL[LL_site][i][lspin_orbit], Model.SpL_On[lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
         }
         for (int i = 0; i < Model.LSpin_LSpin_Jz_site.size() && LL_site - i >= 0; i++) {
            DMRG_Make_Elem_Intersite_Block("LL_LR", Basis.LL, Basis.LR, A_Basis, Inv, System.SzL[LL_site][i][lspin_orbit], Model.SzL_On[lspin_orbit], Model.LSpin_LSpin_Jz_site[i], Ele, "No");
         }
      }
   }
   
   else if (Mat_Type == "LLRL") {
      
      //Interaction: Hopping
      for (int i = 1; i < Model.Ele_Ele_t_site.size() && LL_site - i >= - 1; i++) {
         for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
            DMRG_Make_Elem_Intersite_Block("LL_RL", Basis.LL, Basis.RL, A_Basis, Inv, System.CUp_D  [LL_site][i - 1][ele_orbit], Model.CUp_On    [ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LL_RL", Basis.LL, Basis.RL, A_Basis, Inv, System.CUp    [LL_site][i - 1][ele_orbit], Model.CUp_D_On  [ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LL_RL", Basis.LL, Basis.RL, A_Basis, Inv, System.CDown_D[LL_site][i - 1][ele_orbit], Model.CDown_On  [ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LL_RL", Basis.LL, Basis.RL, A_Basis, Inv, System.CDown  [LL_site][i - 1][ele_orbit], Model.CDown_D_On[ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
         }
      }
      
      //Interaction: Intersite Coulomb
      for (int i = 1; i < Model.Ele_Ele_V_site.size() && LL_site - i >= -1; i++) {
         DMRG_Make_Elem_Intersite_Block("LL_RL", Basis.LL, Basis.RL, A_Basis, Inv, System.NC_Tot[LL_site][i - 1], Model.NC_Tot_On, Model.Ele_Ele_V_site[i], Ele, "No");
      }
      
      //Interaction: Local Spin
      for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
         for (int i = 1; i < Model.LSpin_LSpin_Jxy_site.size() && LL_site - i >= -1; i++) {
            DMRG_Make_Elem_Intersite_Block("LL_RL", Basis.LL, Basis.RL, A_Basis, Inv, System.SpL[LL_site][i - 1][lspin_orbit], Model.SmL_On[lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
            DMRG_Make_Elem_Intersite_Block("LL_RL", Basis.LL, Basis.RL, A_Basis, Inv, System.SmL[LL_site][i - 1][lspin_orbit], Model.SpL_On[lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
         }
         for (int i = 1; i < Model.LSpin_LSpin_Jz_site.size() && LL_site - i >= -1; i++) {
            DMRG_Make_Elem_Intersite_Block("LL_RL", Basis.LL, Basis.RL, A_Basis, Inv, System.SzL[LL_site][i - 1][lspin_orbit], Model.SzL_On[lspin_orbit], Model.LSpin_LSpin_Jz_site[i], Ele, "No");
         }
      }
   }
   
   else if (Mat_Type == "LLRR") {
      
      //Interaction: Hopping
      for (int i = 2; i < Model.Ele_Ele_t_site.size() && LL_site - i >= - 2; i++) {
         for (int j = 0; j < i - 1 && RR_site - i + j >= -2; j++) {
            for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
               DMRG_Make_Elem_Intersite_Block("LL_RR", Basis.LL, Basis.RR, A_Basis, Inv, System.CUp_D  [LL_site][j][ele_orbit], Enviro.CUp    [RR_site][i - j - 2][ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
               DMRG_Make_Elem_Intersite_Block("LL_RR", Basis.LL, Basis.RR, A_Basis, Inv, System.CUp    [LL_site][j][ele_orbit], Enviro.CUp_D  [RR_site][i - j - 2][ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
               DMRG_Make_Elem_Intersite_Block("LL_RR", Basis.LL, Basis.RR, A_Basis, Inv, System.CDown_D[LL_site][j][ele_orbit], Enviro.CDown  [RR_site][i - j - 2][ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
               DMRG_Make_Elem_Intersite_Block("LL_RR", Basis.LL, Basis.RR, A_Basis, Inv, System.CDown  [LL_site][j][ele_orbit], Enviro.CDown_D[RR_site][i - j - 2][ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
            }
         }
      }
      
      //Interaction: Intersite Coulomb
      for (int i = 2; i < Model.Ele_Ele_V_site.size() && LL_site - i >= -2; i++) {
         for (int j = 0; j < i - 1 && RR_site - i + j >= -2; j++) {
            DMRG_Make_Elem_Intersite_Block("LL_RR", Basis.LL, Basis.RR, A_Basis, Inv, System.NC_Tot[LL_site][j], Enviro.NC_Tot[RR_site][i - j - 2], Model.Ele_Ele_V_site[i - j - 2], Ele, "No");
         }
      }
      
      //Interaction: Local Spin
      for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
         for (int i = 2; i < Model.LSpin_LSpin_Jxy_site.size() && LL_site - i >= -2; i++) {
            for (int j = 0; j < i - 1 && RR_site - i + j >= -2; j++) {
               DMRG_Make_Elem_Intersite_Block("LL_RR", Basis.LL, Basis.RR, A_Basis, Inv, System.SpL[LL_site][j][lspin_orbit], Enviro.SmL[RR_site][i - j - 2][lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
               DMRG_Make_Elem_Intersite_Block("LL_RR", Basis.LL, Basis.RR, A_Basis, Inv, System.SmL[LL_site][j][lspin_orbit], Enviro.SpL[RR_site][i - j - 2][lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
            }
         }
         for (int i = 2; i < Model.LSpin_LSpin_Jz_site.size() && LL_site - i >= -2; i++) {
            for (int j = 0; j < i - 1 && RR_site - i + j >= -2; j++) {
               DMRG_Make_Elem_Intersite_Block("LL_RR", Basis.LL, Basis.RR, A_Basis, Inv, System.SzL[LL_site][j][lspin_orbit], Enviro.SzL[RR_site][i - j - 2][lspin_orbit], Model.LSpin_LSpin_Jz_site[i], Ele, "No");
            }
         }
      }
   }
   
   else if (Mat_Type == "LRRL") {
      if (Model.Ele_Ele_t_site.size() > 0) {
         for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
            DMRG_Make_Elem_Intersite_Block("LR_RL", Basis.LR, Basis.RL, A_Basis, Inv, Model.CUp_D_On  [ele_orbit], Model.CUp_On    [ele_orbit], +Model.Ele_Ele_t_site[0], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LR_RL", Basis.LR, Basis.RL, A_Basis, Inv, Model.CUp_On    [ele_orbit], Model.CUp_D_On  [ele_orbit], -Model.Ele_Ele_t_site[0], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LR_RL", Basis.LR, Basis.RL, A_Basis, Inv, Model.CDown_D_On[ele_orbit], Model.CDown_On  [ele_orbit], +Model.Ele_Ele_t_site[0], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LR_RL", Basis.LR, Basis.RL, A_Basis, Inv, Model.CDown_On  [ele_orbit], Model.CDown_D_On[ele_orbit], -Model.Ele_Ele_t_site[0], Ele, "Yes");
         }
      }
      
      if (Model.Ele_Ele_V_site.size() > 0) {
         DMRG_Make_Elem_Intersite_Block("LR_RL", Basis.LR, Basis.RL, A_Basis, Inv, Model.NC_Tot_On, Model.NC_Tot_On, Model.Ele_Ele_V_site[0], Ele, "No");
      }
      
      if (Model.LSpin_LSpin_Jxy_site.size() > 0) {
         for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
            DMRG_Make_Elem_Intersite_Block("LR_RL", Basis.LR, Basis.RL, A_Basis, Inv, Model.SpL_On[lspin_orbit], Model.SmL_On[lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[0], Ele, "No");
            DMRG_Make_Elem_Intersite_Block("LR_RL", Basis.LR, Basis.RL, A_Basis, Inv, Model.SmL_On[lspin_orbit], Model.SpL_On[lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[0], Ele, "No");
         }
      }
      
      if (Model.LSpin_LSpin_Jz_site.size() > 0) {
         for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
            DMRG_Make_Elem_Intersite_Block("LR_RL", Basis.LR, Basis.RL, A_Basis, Inv, Model.SzL_On[lspin_orbit], Model.SzL_On[lspin_orbit], Model.LSpin_LSpin_Jz_site[0], Ele, "No");
         }
      }
   }
   
   else if (Mat_Type == "LRRR") {
      
      //Interaction: Hopping
      for (int i = 1; i < Model.Ele_Ele_t_site.size() && RR_site - i >= - 1; i++) {
         for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
            DMRG_Make_Elem_Intersite_Block("LR_RR", Basis.LR, Basis.RR, A_Basis, Inv, Model.CUp_D_On  [ele_orbit], Enviro.CUp    [RR_site][i - 1][ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LR_RR", Basis.LR, Basis.RR, A_Basis, Inv, Model.CUp_On    [ele_orbit], Enviro.CUp_D  [RR_site][i - 1][ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LR_RR", Basis.LR, Basis.RR, A_Basis, Inv, Model.CDown_D_On[ele_orbit], Enviro.CDown  [RR_site][i - 1][ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("LR_RR", Basis.LR, Basis.RR, A_Basis, Inv, Model.CDown_On  [ele_orbit], Enviro.CDown_D[RR_site][i - 1][ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
         }
      }
      
      //Interaction: Intersite Coulomb
      for (int i = 1; i < Model.Ele_Ele_V_site.size() && LL_site - i >= -1; i++) {
         DMRG_Make_Elem_Intersite_Block("LR_RR", Basis.LR, Basis.RR, A_Basis, Inv, Model.NC_Tot_On, Enviro.NC_Tot[RR_site][i - 1], Model.Ele_Ele_V_site[i], Ele, "No");

      }
      
      //Interaction: Local Spin
      for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
         for (int i = 1; i < Model.LSpin_LSpin_Jxy_site.size() && RR_site - i >= -1; i++) {
            DMRG_Make_Elem_Intersite_Block("LR_RR", Basis.LR, Basis.RR, A_Basis, Inv, Model.SpL_On[lspin_orbit], Enviro.SmL[RR_site][i - 1][lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
            DMRG_Make_Elem_Intersite_Block("LR_RR", Basis.LR, Basis.RR, A_Basis, Inv, Model.SmL_On[lspin_orbit], Enviro.SpL[RR_site][i - 1][lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
         }
         for (int i = 1; i < Model.LSpin_LSpin_Jz_site.size() && LL_site - i >= -1; i++) {
            DMRG_Make_Elem_Intersite_Block("LR_RR", Basis.LR, Basis.RR, A_Basis, Inv, Model.SzL_On[lspin_orbit], Enviro.SzL[RR_site][i - 1][lspin_orbit], Model.LSpin_LSpin_Jz_site[i], Ele, "No");
         }
      }
   }
   
   else if (Mat_Type == "RRRL") {
      
      //Onsite Ham
      DMRG_Make_Elem_Onsite_LLLR_RRRL("RR", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.Ham[RR_site], 1.0);
      DMRG_Make_Elem_Onsite_LLLR_RRRL("RL", Basis.RL, Basis.RR, A_Basis, Inv, Model.Ham_On       , 1.0);
      
      //Interaction: Hopping
      for (int i = 0; i < Model.Ele_Ele_t_site.size() && RR_site - i >= 0; i++) {
         for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
            DMRG_Make_Elem_Intersite_Block("RR_RL", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.CUp_D  [RR_site][i][ele_orbit], Model.CUp_On    [ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("RR_RL", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.CUp    [RR_site][i][ele_orbit], Model.CUp_D_On  [ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("RR_RL", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.CDown_D[RR_site][i][ele_orbit], Model.CDown_On  [ele_orbit], +Model.Ele_Ele_t_site[i], Ele, "Yes");
            DMRG_Make_Elem_Intersite_Block("RR_RL", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.CDown  [RR_site][i][ele_orbit], Model.CDown_D_On[ele_orbit], -Model.Ele_Ele_t_site[i], Ele, "Yes");
         }
      }
      
      //Interaction: Intersite Coulomb
      for (int i = 0; i < Model.Ele_Ele_V_site.size() && RR_site - i >= 0; i++) {
         DMRG_Make_Elem_Intersite_Block("RR_RL", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.NC_Tot[RR_site][i], Model.NC_Tot_On, Model.Ele_Ele_V_site[i], Ele, "No");
      }
      
      //Interaction: Local Spin
      for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
         for (int i = 0; i < Model.LSpin_LSpin_Jxy_site.size() && RR_site - i >= 0; i++) {
            DMRG_Make_Elem_Intersite_Block("RR_RL", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.SpL[RR_site][i][lspin_orbit], Model.SmL_On[lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
            DMRG_Make_Elem_Intersite_Block("RR_RL", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.SmL[RR_site][i][lspin_orbit], Model.SpL_On[lspin_orbit], 0.5*Model.LSpin_LSpin_Jxy_site[i], Ele, "No");
         }
         for (int i = 0; i < Model.LSpin_LSpin_Jz_site.size() && LL_site - i >= 0; i++) {
            DMRG_Make_Elem_Intersite_Block("RR_RL", Basis.RR, Basis.RL, A_Basis, Inv, Enviro.SzL[RR_site][i][lspin_orbit], Model.SzL_On[lspin_orbit], Model.LSpin_LSpin_Jz_site[i], Ele, "No");
         }
      }
   }
   
   else {
      std::cout << "Error in Make_Elem_Ham" << std::endl;
      std::exit(0);
   }
   
   DMRG_Clear_Check_Basis(A_Basis, Inv);
   
}
