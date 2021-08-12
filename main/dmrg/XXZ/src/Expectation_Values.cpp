//
//  Expectation_Values.cpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/06/30.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Expectation_Values(DMRG_Ground_State &GS, DMRG_Basis_LLLRRRRL &Bases_LLLRRRRL, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_XXZ &Model, DMRG_Param &Dmrg_Param) {
   
   int LL_site = Block.LL_site;
   int RR_site = Block.RR_site;
   
   int c1 = (Dmrg_Param.now_sweep == 0);
   int c2 = (Dmrg_Param.tot_sweep - Dmrg_Param.now_sweep >= 2);
   int c3 = (LL_site + RR_site + 4 != Model.system_size);
   int c4 = (LL_site != RR_site);
   
   if (c1 || c2 || c3 || c4) {
      return;
   }
   
   //Output_Average_Values(GS.val, "avg_energy.txt", LL_site, RR_site, Dmrg_Param, Model);

   std::vector<CRS> Mat_System, Mat_Enviro, Mat_CF;
   std::vector<double> Sz, SxSx, SySy, SzSz;
   std::vector<double> Sx_CF, Sy_CF, Sz_CF;
   std::vector<int> Dummy;

   //Sz
   DMRG_Transform_Matrix_One(Model.Sz_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.Sz_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(Sz, Mat_System, Model.Sz_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Model.Output_Onsite_Values(Sz, "Sz.txt");
   
   //Sz_CF
   DMRG_Transform_Matrix_Two(Model.Sz_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(Sz_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.Sz_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Model.Output_Intersite_Values(Sz_CF, "Sz_CF.txt");
   
   //SxSx
   DMRG_Transform_Matrix_One(Model.SxSx_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxSx_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxSx, Mat_System, Model.SxSx_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Model.Output_Onsite_Values(SxSx, "SxSx.txt");
   
   //SySy
   DMRG_Transform_Matrix_One(Model.SySy_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SySy_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SySy, Mat_System, Model.SySy_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Model.Output_Onsite_Values(SySy, "SySy.txt");
   
   //SzSz
   DMRG_Transform_Matrix_One(Model.SzSz_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzSz_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzSz, Mat_System, Model.SzSz_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Model.Output_Onsite_Values(SzSz, "SzSz.txt");
   
   std::vector<DMRG_Basis_LLLRRRRL> W_Basis(2);
   W_Basis[0].qn1 = Model.tot_sz + 2;
   DMRG_Get_Basis_LLLRRRRL(W_Basis[0], Basis_System, Basis_Enviro, Block, Model.num_of_qn, "No", Model.p_thread);
   W_Basis[1].qn1 = Model.tot_sz - 2;
   DMRG_Get_Basis_LLLRRRRL(W_Basis[1], Basis_System, Basis_Enviro, Block, Model.num_of_qn, "No", Model.p_thread);
   
   //Sx_CF
   DMRG_Transform_Matrix_One(Model.Sx_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.Sx_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.Sx_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(Sx_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.Sx_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Model.Output_Intersite_Values(Sx_CF, "Sx_CF.txt");
   
   //Sy_CF
   DMRG_Transform_Matrix_One(Model.iSy_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.iSy_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.iSy_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(Sy_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.iSy_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Model.Output_Intersite_Values(Sy_CF, "Sy_CF.txt");
   
}
