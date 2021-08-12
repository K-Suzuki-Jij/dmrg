
//
//  Expectation_Values.cpp
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/07/10.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Expectation_Values(DMRG_Ground_State &GS, DMRG_Basis_LLLRRRRL &Bases_LLLRRRRL, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_AKLM &Model, DMRG_Param &Dmrg_Param) {
   
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
   std::vector<double> SCSL, Ham;
   std::vector<double> SzL, SxLSxL, SyLSyL, SzLSzL;
   std::vector<double> SzC, SxCSxC, SyCSyC, SzCSzC;
   std::vector<double> NCUp, NCDown, NC;
   std::vector<double> SxL_CF, SyL_CF, SzL_CF;
   std::vector<double> SxC_CF, SyC_CF, SzC_CF;
   std::vector<double> NCUp_CF, NCDown_CF, NC_CF;
   std::vector<int> Dummy;
   
   DMRG_Transform_Matrix_One(Model.SCSL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SCSL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SCSL, Mat_System, Model.SCSL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SCSL, "SCSL.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.Ham_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.Ham_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(Ham, Mat_System, Model.Ham_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(Ham, "Ham.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzL, Mat_System, Model.SzL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzL, "SzL.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.SzL_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SzL_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SzL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SzL_CF, SzL, "SzL_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SxLSxL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxLSxL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxLSxL, Mat_System, Model.SxLSxL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SxLSxL, "SxLSxL.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SyLSyL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SyLSyL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SyLSyL, Mat_System, Model.SyLSyL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SyLSyL, "SyLSyL.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzLSzL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzLSzL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzLSzL, Mat_System, Model.SzLSzL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzLSzL, "SzLSzL.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzC, Mat_System, Model.SzC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzC, "SzC.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.SzC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SzC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SzC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SzC_CF, SzC, "SzC_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SxCSxC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxCSxC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxCSxC, Mat_System, Model.SxCSxC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SxCSxC, "SxCSxC.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SyCSyC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SyCSyC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SyCSyC, Mat_System, Model.SyCSyC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SyCSyC, "SyCSyC.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzCSzC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzCSzC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzCSzC, Mat_System, Model.SzCSzC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzCSzC, "SzCSzC.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NCUp_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCUp_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCUp, Mat_System, Model.NCUp_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NCUp, "NCUp.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NCUp_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCUp_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCUp_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NCUp_CF, NCUp, "NCUp_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NCDown_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCDown_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCDown, Mat_System, Model.NCDown_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NCDown, "NCDown.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NCDown_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCDown_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCDown_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NCDown_CF, NCDown, "NCDown_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NC, Mat_System, Model.NC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NC, "NC.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NC_CF, NC, "NC_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   std::vector<DMRG_Basis_LLLRRRRL> W_Basis(2);
   W_Basis[0].qn1 = Model.tot_ele;
   W_Basis[0].qn2 = Model.tot_sz + 2;
   DMRG_Get_Basis_LLLRRRRL(W_Basis[0], Basis_System, Basis_Enviro, Block, Model.num_of_qn, "No", Model.p_thread);
   W_Basis[1].qn1 = Model.tot_ele;
   W_Basis[1].qn2 = Model.tot_sz - 2;
   DMRG_Get_Basis_LLLRRRRL(W_Basis[1], Basis_System, Basis_Enviro, Block, Model.num_of_qn, "No", Model.p_thread);
   
   std::vector<double> SxL(Model.system_size, 0.0), SyL(Model.system_size, 0.0);
   std::vector<double> SxC(Model.system_size, 0.0), SyC(Model.system_size, 0.0);
   
   DMRG_Transform_Matrix_One(Model.SxL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.SxL_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SxL_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SxL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SxL_CF, SxL, "SxL_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.iSyL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.iSyL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.iSyL_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SyL_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.iSyL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SyL_CF, SyL, "SyL_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SxC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.SxC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SxC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SxC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SxC_CF, SxC, "SxC_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.iSyC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.iSyC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.iSyC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SyC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.iSyC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SyC_CF, SyC, "SyC_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
}
