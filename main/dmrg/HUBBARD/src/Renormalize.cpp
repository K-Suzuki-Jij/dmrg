//
//  Renormalize.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Renormalize(Block_Operator &System, Block_Operator &Enviro, DMRG_Basic_Information &Info, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, Model_1D_HUBBARD &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param) {
   
   double start = omp_get_wtime();
   Dmrg_Param.renorm_now_iter++;
   
   DMRG_Time Time;
   
   if (Info.Block_Compo.LL_site == 0) {
      Get_Operator_Edge(System, Info.Block_Compo, Basis_System, Model);
   }
   
   Find_GS_QN(Info.QN, Info.Block_Compo, Model);
   
   DMRG_Construct_Superblock(Info, Basis_System, Basis_Enviro, Model.num_of_qn, Dmrg_Param.Initial_Guess, "No", Model.p_thread, Time);
   
   DMRG_Initial_Guess(Info.GS.Vector, Info.GS.Vector_Guess, Info.Block_Compo, Info.Basis_LLLRRRRL, Info.Basis_LLLR, Basis_System, Basis_Enviro, "Ele", Dmrg_Param.Initial_Guess ,Model.p_thread);

   Diagonalize_Ham_LLLRRRRL(System, Enviro, Info, Basis_System, Basis_Enviro, Model, Dmrg_Param, Diag_Param, Time);
   
   DMRG_T_Mat T_Mat;
   DMRG_Get_Trans_Matrix(T_Mat, Info, Model.num_of_qn, Dmrg_Param.max_dim_system, Model.p_thread, Time);
   
   Renormalize_System(System, Info.Basis_LLLR, Info.Block_Compo, Basis_System, T_Mat, Model);
   
   DMRG_Store_Trans_Matrix(T_Mat.Trans_Mat, Basis_System, Info.Basis_LLLR, Dmrg_Param, Info.Block_Compo.LL_site, Model.system_size);
   
   Expectation_Values(Info.GS.Vector, Info.Basis_LLLRRRRL, Info.Block_Compo, Basis_System, Basis_Enviro, Info.QN, Model, Dmrg_Param);
   
   Time.total = omp_get_wtime() - start;
   
   DMRG_Print_Status("Ele+Spin", Model.BC, Info, Dmrg_Param, Time);
   DMRG_Output_Time(Time);
   
}

