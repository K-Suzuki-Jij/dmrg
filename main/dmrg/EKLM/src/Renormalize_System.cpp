//
//  Renormalize_System.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/08.
//

#include "Header.hpp"

void Renormalize_System(Block_Operator &System, const DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, const DMRG_Block_Information &Block, const DMRG_T_Mat &T_Mat, const Model_1D_EKLM &Model) {
   
   int LL_site = Block.LL_site;
   int dim_LL  = Block.dim_LL;
   std::vector<int> Ele_LL(dim_LL, 0);
   
#pragma omp parallel for num_threads (Model.p_threads)
   for (int i = 0; i < dim_LL; i++) {
      for (int ele_row = Basis_System.qn_LL_LL_stored_ele_start; ele_row < Basis_System.qn_LL_LL_stored_ele_end; ele_row++) {
         Ele_LL[i] += Basis_System.QN_LL_LL_Stored[LL_site][ele_row][i];
      }
   }
      
   CRS Ham_LLLR;
   Get_Ham_Block("LLLR", System, Block_Operator(), Basis, Ele_LL, Block, Model, Ham_LLLR);
   Change_Matrix_Basis(Ham_LLLR, T_Mat.Trans_Mat, System.Ham[LL_site + 1], Model.p_threads);

   
   //LR Renormalize
   for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
#pragma omp parallel sections num_threads (Model.p_threads)
      {
#pragma omp section
      DMRG_Renormalize_Matrix("LR", Model.CUp_On  [ele_orbit], T_Mat.Trans_Mat, System.CUp  [LL_site + 1][0][ele_orbit], Block, "Yes", Basis.LLLR, Ele_LL, 1);
#pragma omp section
      DMRG_Renormalize_Matrix("LR", Model.CDown_On[ele_orbit], T_Mat.Trans_Mat, System.CDown[LL_site + 1][0][ele_orbit], Block, "Yes", Basis.LLLR, Ele_LL, 1);
      }
   }
   
   DMRG_Renormalize_Matrix("LR", Model.NC_Tot_On, T_Mat.Trans_Mat, System.NC_Tot[LL_site + 1][0], Block, "No", Basis.LLLR, std::vector<int>(), Model.p_threads);
   
   for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
#pragma omp parallel sections num_threads (Model.p_threads)
      {
#pragma omp section
      DMRG_Renormalize_Matrix("LR", Model.SpL_On[lspin_orbit], T_Mat.Trans_Mat, System.SpL[LL_site + 1][0][lspin_orbit], Block, "No", Basis.LLLR, std::vector<int>(), 1);
#pragma omp section
      DMRG_Renormalize_Matrix("LR", Model.SzL_On[lspin_orbit], T_Mat.Trans_Mat, System.SzL[LL_site + 1][0][lspin_orbit], Block, "No", Basis.LLLR, std::vector<int>(), 1);
      }
   }

   for (int ele_orbit = 0; ele_orbit < Model.num_ele_orbit; ele_orbit++) {
#pragma omp parallel sections num_threads (Model.p_threads)
      {
#pragma omp section
      Matrix_Transpose(System.CUp  [LL_site + 1][0][ele_orbit], System.CUp_D  [LL_site + 1][0][ele_orbit]);
#pragma omp section
      Matrix_Transpose(System.CDown[LL_site + 1][0][ele_orbit], System.CDown_D[LL_site + 1][0][ele_orbit]);
      }
   }
   
   for (int lspin_orbit = 0; lspin_orbit < Model.num_lspin_orbit; lspin_orbit++) {
      Matrix_Transpose(System.SpL[LL_site + 1][0][lspin_orbit], System.SmL[LL_site + 1][0][lspin_orbit]);
   }
   
   
   int dim_renorm = T_Mat.Trans_Mat.col_dim;
   
   Basis_System.QN_LL_LL_Stored[LL_site + 1].resize(Model.num_of_qn);
   for (int i = 0; i < Model.num_of_qn; i++) {
      Basis_System.QN_LL_LL_Stored[LL_site + 1][i].resize(dim_renorm);
      for (int j = 0; j < dim_renorm; j++) {
         Basis_System.QN_LL_LL_Stored[LL_site + 1][i][j] = T_Mat.QN[i][j];
      }
   }
   
}
