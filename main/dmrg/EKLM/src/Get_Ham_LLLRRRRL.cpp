//
//  Get_Ham_LLLRRRRL.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/04.
//

#include "Header.hpp"

void Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, const Block_Operator &System, const Block_Operator &Enviro, const DMRG_Basis &Basis, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, const DMRG_Block_Information &Block, DMRG_Time &Time, const Model_1D_EKLM &Model) {
   
   double start = omp_get_wtime();
   
   int LL_site           = Block.LL_site;
   int RR_site           = Block.RR_site;
   int dim_onsite        = Block.dim_onsite;
   int dim_LL            = Block.dim_LL;
   int dim_RR            = Block.dim_RR;
   
   DMRG_Block_Hamiltonian Block_Ham;

   Block_Ham.Ele_LL.resize(dim_LL);
   Block_Ham.Ele_RR.resize(dim_RR);
   Block_Ham.Ele_On.resize(dim_onsite);
   
#pragma omp parallel for num_threads (Model.p_threads)
   for (int i = 0; i < dim_LL; i++) {
      Block_Ham.Ele_LL[i] = 0;
   }
#pragma omp parallel for num_threads (Model.p_threads)
   for (int i = 0; i < dim_RR; i++) {
      Block_Ham.Ele_RR[i] = 0;
   }
#pragma omp parallel for num_threads (Model.p_threads)
   for (int i = 0; i < dim_onsite; i++) {
      Block_Ham.Ele_On[i] = 0;
   }
   
#pragma omp parallel for num_threads (Model.p_threads)
   for (int i = 0; i < dim_LL; i++) {
      for (int ele_row = Basis_System.qn_LL_LL_stored_ele_start; ele_row < Basis_System.qn_LL_LL_stored_ele_end; ele_row++) {
         Block_Ham.Ele_LL[i] += Basis_System.QN_LL_LL_Stored[LL_site][ele_row][i];
      }
   }
   
#pragma omp parallel for num_threads (Model.p_threads)
   for (int i = 0; i < dim_RR; i++) {
      for (int ele_row = Basis_Enviro.qn_LL_LL_stored_ele_start; ele_row < Basis_Enviro.qn_LL_LL_stored_ele_end; ele_row++) {
         Block_Ham.Ele_RR[i] += Basis_Enviro.QN_LL_LL_Stored[RR_site][ele_row][i];
      }
   }
   
#pragma omp parallel for num_threads (Model.p_threads)
   for (int i = 0; i < dim_onsite; i++) {
      for (int ele_row = Basis_System.qn_LL_LL_stored_ele_start; ele_row < Basis_System.qn_LL_LL_stored_ele_end; ele_row++) {
         Block_Ham.Ele_On[i] += Basis_System.QN_LL_LL_Stored[0][ele_row][i];
      }
   }

   Get_Ham_Block("LLLR", System          , Block_Operator(), Basis, Block_Ham.Ele_LL, Block, Model, Block_Ham.LLLR);
   Get_Ham_Block("LLRR", System          , Enviro          , Basis, Block_Ham.Ele_LL, Block, Model, Block_Ham.LLRR);
   Get_Ham_Block("LLRL", System          , Block_Operator(), Basis, Block_Ham.Ele_LL, Block, Model, Block_Ham.LLRL);
   Get_Ham_Block("LRRR", Block_Operator(), Enviro          , Basis, Block_Ham.Ele_On, Block, Model, Block_Ham.LRRR);
   Get_Ham_Block("LRRL", Block_Operator(), Block_Operator(), Basis, Block_Ham.Ele_On, Block, Model, Block_Ham.LRRL);
   Get_Ham_Block("RRRL", Block_Operator(), Enviro          , Basis, Block_Ham.Ele_RR, Block, Model, Block_Ham.RRRL);
   
   std::string Sign_Flag = "No";
   
   if (Basis_System.qn_LL_LL_stored_ele_start < Basis_System.qn_LL_LL_stored_ele_end) {
      Sign_Flag = "Yes";
   }
   
   DMRG_Get_Ham_LLLRRRRL(Ham_LLLRRRRL, Block_Ham, Basis, Block, Sign_Flag, Model.p_threads);

   Time.make_ham = omp_get_wtime() - start;
   
}
