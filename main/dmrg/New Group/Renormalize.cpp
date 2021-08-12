//
//  Renormalize.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/04.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Renormalize(Block_Operator &System, Block_Operator &Enviro, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, DMRG_Basis &Basis, DMRG_Ground_State &GS, Model_1D_TAKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param) {
   
   double start = omp_get_wtime();
   Dmrg_Param.renorm_now_iter++;
   
   DMRG_Time Time;
   DMRG_T_Mat T_Mat;

   Find_GS_QN(Basis.LLLRRRRL.qn1, Basis.LLLRRRRL.qn2, Basis.LLLRRRRL.qn3, Block, Model);

   DMRG_Construct_Superblock(Basis, Basis_System, Basis_Enviro, Block, Model.num_of_qn, Dmrg_Param.Initial_Guess, "No", Model.p_thread, Time);
      
   DMRG_Initial_Guess(GS, Basis.LLLRRRRL, Basis.LLLR, Basis_System, Basis_Enviro, Block, "Ele+Ele", Dmrg_Param.Initial_Guess ,Model.p_thread);
   
   CRS Ham_LLLRRRRL;
   Get_Ham_LLLRRRRL(Ham_LLLRRRRL, System, Enviro, Basis, Basis_System, Basis_Enviro, Block, Time, Model);
   DMRG_Diagonalize_Ham_LLLRRRRL(Ham_LLLRRRRL, GS, Basis, Block, Dmrg_Param, Diag_Param, Time, Model.p_thread);
   Free_CRS(Ham_LLLRRRRL);
   
   DMRG_Get_Trans_Matrix(T_Mat, Basis.LLLRRRRL, Basis.LLLR, GS, Model.num_of_qn, Dmrg_Param.max_dim_system, Model.p_thread, Time);
   
   Renormalize_System(System, Basis, Basis_System, Block.LL_site, T_Mat, Model);
   
   DMRG_Store_Trans_Matrix(T_Mat.Trans_Mat, Basis_System, Basis.LLLR, Dmrg_Param, Block.LL_site, Model.system_size);
   
   Expectation_Values(GS, Basis.LLLRRRRL, Basis_System, Basis_Enviro, Block, Model, Dmrg_Param);
      
   Time.total = omp_get_wtime() - start;
   
   DMRG_Print_Status("Ele+Ele+Spin", Model.BC, GS, Basis.LLLRRRRL, Basis_System, Basis_Enviro, Block, Dmrg_Param, Time);

   DMRG_Output_Time(Time);
   
   
}

