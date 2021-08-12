//
//  Expectation_Values.cpp
//  1D_AKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2020/08/01.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Expectation_Values(DMRG_Ground_State &GS, DMRG_LLLRRRRL_Basis &Bases_LLLRRRRL, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Quntum_Number &QN, Model_1D_AKLM_TVF &Model, DMRG_Param &Dmrg_Param) {
   
   int c1 = (Dmrg_Param.now_sweep == 0);
   int c2 = (Dmrg_Param.tot_sweep - Dmrg_Param.now_sweep >= 2);
   int c3 = (Block_Compo.LL_site + Block_Compo.RR_site + 4 != Model.system_size);
   int c4 = (Block_Compo.LL_site != Block_Compo.RR_site);
   
   if (c1 || c2 || c3 || c4) {
      return;
   }
   
   Output_Average_Values(GS.val, "avg_energy.txt", Block_Compo, Dmrg_Param, Model);
   
   int LL_site = Block_Compo.LL_site;
   int RR_site = Block_Compo.RR_site;
   
   std::vector<CRS> Mat_System, Mat_Enviro, Mat_CF;
   std::vector<double> SCSL, Ham;
   std::vector<double> SxL, SxLSxL, SyLSyL, SzLSzL;
   std::vector<double> SxC, SxCSxC, SyCSyC, SzCSzC;
   std::vector<double> NCEven, NCOdd, NC;
   std::vector<double> SxL_CF, SyL_CF, SzL_CF;
   std::vector<double> SxC_CF, SyC_CF, SzC_CF;
   std::vector<double> NCEven_CF, NCOdd_CF, NC_CF;
   std::vector<int> Dummy;
   
   DMRG_Transform_Matrix_One(Model.SCSL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SCSL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SCSL, Mat_System, Model.SCSL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(SCSL     , "SCSL.txt"    , Block_Compo, Dmrg_Param, Model);
   Output_Average_Values(SCSL    , "avg_SCSL.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SCSL, "Fourier_SCSL.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SxL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxL, Mat_System, Model.SxL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(SxL, "SxL.txt", Block_Compo, Dmrg_Param, Model);
   Output_Average_Values(SxL, "avg_SxL.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SxL, "Fourier_SxL.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_Two(Model.SxL_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SxL_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SxL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(SxL_CF, SxL, "SxL_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SxL_CF, "Fourier_SxL_CF.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.SxLSxL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxLSxL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxLSxL, Mat_System, Model.SxLSxL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(SxLSxL, "SxLSxL.txt", Block_Compo, Dmrg_Param, Model);
   Output_Average_Values(SxLSxL, "avg_SxLSxL.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SxLSxL, "Fourier_SxLSxL.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.SxC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxC, Mat_System, Model.SxC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(SxC, "SxC.txt", Block_Compo, Dmrg_Param, Model);
   Output_Average_Values(SxC, "avg_SxC.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SxC, "Fourier_SxC.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_Two(Model.SxC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SxC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SxC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(SxC_CF, SxC, "SxC_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SxC_CF, "Fourier_SxC_CF.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.SxCSxC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxCSxC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxCSxC, Mat_System, Model.SxCSxC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(SxCSxC, "SxCSxC.txt", Block_Compo, Dmrg_Param, Model);
   Output_Average_Values(SxCSxC, "avg_SxCSxC.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SxCSxC, "Fourier_SxCSxC.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.NCEven_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCEven_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCEven, Mat_System, Model.NCEven_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(NCEven, "NCEven.txt", Block_Compo, Dmrg_Param, Model);
   Output_Average_Values(NCEven, "avg_NCEven.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(NCEven, "Fourier_NCEven.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_Two(Model.NCEven_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCEven_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCEven_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(NCEven_CF, NCEven, "NCEven_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(NCEven_CF, "Fourier_NCEven_CF.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.NCOdd_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCOdd_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCOdd, Mat_System, Model.NCOdd_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(NCOdd, "NCOdd.txt", Block_Compo, Dmrg_Param, Model);
   Output_Average_Values(NCOdd, "avg_NCOdd.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(NCOdd, "Fourier_NCOdd.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_Two(Model.NCOdd_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCOdd_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCOdd_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(NCOdd_CF, NCOdd, "NCOdd_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(NCOdd_CF, "Fourier_NCOdd_CF.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.NC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NC, Mat_System, Model.NC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(NC, "NC.txt", Block_Compo, Dmrg_Param, Model);
   Output_Average_Values(NC, "avg_NC.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(NC, "Fourier_NC.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_Two(Model.NC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(NC_CF, NC, "NC_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(NC_CF, "Fourier_NC_CF.txt", Block_Compo, Dmrg_Param, Model);

   std::vector<DMRG_LLLRRRRL_Basis> W_Basis(2);
   DMRG_Quntum_Number W_QN;
   W_QN.qn1_LLLRRRRL = Model.tot_ele;
   W_QN.qn2_LLLRRRRL = Model.tot_parity;
   DMRG_Get_Basis_LLLRRRRL(W_QN, W_Basis[0], Block_Compo, Basis_System, Basis_Enviro, Model.num_of_qn, "qn2", Model.p_thread);
   W_QN.qn1_LLLRRRRL = Model.tot_ele;
   W_QN.qn2_LLLRRRRL = (Model.tot_parity + 1)%2;
   DMRG_Get_Basis_LLLRRRRL(W_QN, W_Basis[1], Block_Compo, Basis_System, Basis_Enviro, Model.num_of_qn, "qn2", Model.p_thread);
   
   std::vector<double> SzL(Model.system_size, 0.0), SyL(Model.system_size, 0.0), SzC(Model.system_size, 0.0), SyC(Model.system_size, 0.0);
   
   DMRG_Transform_Matrix_One(Model.SzL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.SzL_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SzL_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SzL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, Block_Compo, Model.p_thread);
   Output_Intersite_Values(SzL_CF, SzL, "SzL_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SzL_CF, "Fourier_SzL_CF.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.iSyL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.iSyL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.iSyL_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SyL_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.iSyL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, Block_Compo, Model.p_thread);
   Output_Intersite_Values(SyL_CF, SyL, "SyL_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SyL_CF, "Fourier_SyL_CF.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.SzC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.SzC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SzC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SzC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, Block_Compo, Model.p_thread);
   Output_Intersite_Values(SzC_CF, SzC, "SzC_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SzC_CF, "Fourier_SzC_CF.txt", Block_Compo, Dmrg_Param, Model);

   DMRG_Transform_Matrix_One(Model.iSyC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.iSyC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.iSyC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SyC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.iSyC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, Block_Compo, Model.p_thread);
   Output_Intersite_Values(SyC_CF, SyC, "SyC_CF.txt", Block_Compo, Dmrg_Param, Model);
   Output_Fourier_Components(SyC_CF, "Fourier_SyC_CF.txt", Block_Compo, Dmrg_Param, Model);

}
