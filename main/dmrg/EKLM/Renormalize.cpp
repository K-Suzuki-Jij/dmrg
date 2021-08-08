//
//  Renormalize.cpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/01.
//

#include "Header.hpp"

void Renormalize(Block_Operator &System, Block_Operator &Enviro, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, DMRG_Basis &Basis, DMRG_Ground_State &GS, Model_1D_EKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param) {
   
   double start = omp_get_wtime();
   Dmrg_Param.renorm_now_iter++;
   
   DMRG_Time Time;
   DMRG_T_Mat T_Mat;
   
   Find_GS_QN(Basis.LLLRRRRL, Block, Model);
   
   DMRG_Construct_Superblock(Basis, Basis_System, Basis_Enviro, Block, Dmrg_Param.Initial_Guess, Time, Model.p_threads);
   
   DMRG_Initial_Guess(GS, Basis.LLLRRRRL, Basis.LLLR, Basis_System, Basis_Enviro, Block, Dmrg_Param.Initial_Guess, Model.p_threads);
   
   CRS Ham_LLLRRRRL;
   Get_Ham_LLLRRRRL(Ham_LLLRRRRL, System, Enviro, Basis, Basis_System, Basis_Enviro, Block, Time, Model);
   
   DMRG_Diagonalize_Ham_LLLRRRRL(Ham_LLLRRRRL, GS, Basis, Block, Dmrg_Param, Diag_Param, Time, Model.p_threads);
   Free_CRS(Ham_LLLRRRRL);
      
   DMRG_Get_Trans_Matrix(T_Mat, Basis.LLLRRRRL, Basis.LLLR, GS, Model.num_of_qn, Dmrg_Param.max_dim_system, Model.p_threads, Time);

   Renormalize_System(System, Basis, Basis_System, Block, T_Mat, Model);
   
   DMRG_Store_Trans_Matrix(T_Mat.Trans_Mat, Basis_System, Basis.LLLR, Dmrg_Param, Block.LL_site, Model.system_size);

   //Expectation_Values(GS, Basis.LLLRRRRL, Basis_System, Basis_Enviro, Block, Model, Dmrg_Param);
   
   Time.total = omp_get_wtime() - start;

   DMRG_Print_Status(Model.num_of_qn, Model.BC, GS, Basis.LLLRRRRL, Basis_System, Basis_Enviro, Block, Dmrg_Param, Time);
   
}
