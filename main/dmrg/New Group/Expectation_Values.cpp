//
//  Expectation_Values.cpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/11.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Expectation_Values(DMRG_Ground_State &GS, DMRG_Basis_LLLRRRRL &Bases_LLLRRRRL, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_TAKLM &Model, DMRG_Param &Dmrg_Param) {
   
   int LL_site = Block.LL_site;
   int RR_site = Block.RR_site;
   
   int c1 = (Dmrg_Param.now_sweep == 0);
   int c2 = (Dmrg_Param.tot_sweep - Dmrg_Param.now_sweep >= 2);
   int c3 = (LL_site + RR_site + 4 != Model.system_size);
   int c4 = (LL_site != RR_site);
   
   if (c1 || c2 || c3 || c4) {
      return;
   }
   
   Output_Average_Values(GS.val, "avg_energy.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   std::vector<CRS> Mat_System, Mat_Enviro, Mat_CF;
   std::vector<double> SC_1SL, SC_2SL;
   std::vector<double> SzL, SxLSxL, SyLSyL, SzLSzL;
   std::vector<double> SzC_1, SxC_1SxC_1, SyC_1SyC_1, SzC_1SzC_1;
   std::vector<double> SzC_2, SxC_2SxC_2, SyC_2SyC_2, SzC_2SzC_2;
   std::vector<double> NCUp_1, NCDown_1, NC_1;
   std::vector<double> NCUp_2, NCDown_2, NC_2, NC;
   std::vector<double> SxL_CF, SyL_CF, SzL_CF;
   std::vector<double> SxC_1_CF, SyC_1_CF, SzC_1_CF;
   std::vector<double> SxC_2_CF, SyC_2_CF, SzC_2_CF;
   std::vector<double> NCUp_1_CF, NCDown_1_CF, NC_1_CF;
   std::vector<double> NCUp_2_CF, NCDown_2_CF, NC_2_CF, NC_CF;
   std::vector<int> Dummy;
   
   double s = omp_get_wtime();
   DMRG_Transform_Matrix_One(Model.SC_1SL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   printf("SCSL1=%lf\n",omp_get_wtime() - s);
   s = omp_get_wtime();
   DMRG_Transform_Matrix_One(Model.SC_1SL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   printf("SCSL2=%lf\n",omp_get_wtime() - s);
   s = omp_get_wtime();
   DMRG_Expectation_Onsite(SC_1SL, Mat_System, Model.SC_1SL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   printf("SCSL3=%lf\n",omp_get_wtime() - s);
   Output_Onsite_Values(SC_1SL, "SC_1SL.txt", LL_site, RR_site, Dmrg_Param, Model);
   Output_Average_Values(SC_1SL, "avg_SC_1SL.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SC_2SL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SC_2SL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SC_2SL, Mat_System, Model.SC_2SL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SC_2SL, "SC_2SL.txt", LL_site, RR_site, Dmrg_Param, Model);
   Output_Average_Values(SC_2SL, "avg_SC_2SL.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzL_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzL_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzL, Mat_System, Model.SzL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzL, "SzL.txt", LL_site, RR_site, Dmrg_Param, Model);
   Output_Average_Values(SzL, "avg_SzL.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   s = omp_get_wtime();
   DMRG_Transform_Matrix_Two(Model.SzL_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   printf("SzLCF1=%lf\n",omp_get_wtime() - s);
   s = omp_get_wtime();
   DMRG_Expectation_Intersite(SzL_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SzL_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   printf("SzLCF2=%lf\n",omp_get_wtime() - s);
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
   
   DMRG_Transform_Matrix_One(Model.SzC_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzC_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzC_1, Mat_System, Model.SzC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzC_1, "SzC_1.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.SzC_1_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SzC_1_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SzC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SzC_1_CF, SzC_1, "SzC_1_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SxC_1SxC_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxC_1SxC_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxC_1SxC_1, Mat_System, Model.SxC_1SxC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SxC_1SxC_1, "SxC_1SxC_1.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SyC_1SyC_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SyC_1SyC_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SyC_1SyC_1, Mat_System, Model.SyC_1SyC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SyC_1SyC_1, "SyC_1SyC_1.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzC_1SzC_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzC_1SzC_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzC_1SzC_1, Mat_System, Model.SzC_1SzC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzC_1SzC_1, "SzC_1SzC_1.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzC_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzC_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzC_2, Mat_System, Model.SzC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzC_2, "SzC_2.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.SzC_2_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SzC_2_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SzC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SzC_2_CF, SzC_2, "SzC_2_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SxC_2SxC_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxC_2SxC_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxC_2SxC_2, Mat_System, Model.SxC_2SxC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SxC_2SxC_2, "SxC_2SxC_2.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SyC_2SyC_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SyC_2SyC_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SyC_2SyC_2, Mat_System, Model.SyC_2SyC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SyC_2SyC_2, "SyC_2SyC_2.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzC_2SzC_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzC_2SzC_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzC_2SzC_2, Mat_System, Model.SzC_2SzC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(SzC_2SzC_2, "SzC_2SzC_2.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NCUp_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCUp_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCUp_1, Mat_System, Model.NCUp_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NCUp_1, "NCUp_1.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NCUp_1_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCUp_1_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCUp_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NCUp_1_CF, NCUp_1, "NCUp_1_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NCDown_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCDown_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCDown_1, Mat_System, Model.NCDown_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NCDown_1, "NCDown_1.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NCDown_1_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCDown_1_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCDown_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NCDown_1_CF, NCDown_1, "NCDown_1_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NC_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NC_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NC_1, Mat_System, Model.NC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NC_1, "NC_1.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NC_1_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NC_1_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NC_1_CF, NC_1, "NC_1_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NCUp_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCUp_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCUp_2, Mat_System, Model.NCUp_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NCUp_2, "NCUp_2.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NCUp_2_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCUp_2_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCUp_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NCUp_2_CF, NCUp_2, "NCUp_2_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NCDown_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCDown_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCDown_2, Mat_System, Model.NCDown_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NCDown_2, "NCDown_2.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NCDown_2_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCDown_2_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCDown_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NCDown_2_CF, NCDown_2, "NCDown_2_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NC_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NC_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NC_2, Mat_System, Model.NC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NC_2, "NC_2.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NC_2_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NC_2_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NC_2_CF, NC_2, "NC_2_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NC, Mat_System, Model.NC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Onsite_Values(NC, "NC.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NC_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(NC_CF, NC, "NC_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   std::vector<DMRG_Basis_LLLRRRRL> W_Basis(2);
   W_Basis[0].qn1 = Model.tot_ele_1;
   W_Basis[0].qn2 = Model.tot_ele_2;
   W_Basis[0].qn3 = Model.tot_sz + 2;
   DMRG_Get_Basis_LLLRRRRL(W_Basis[0], Basis_System, Basis_Enviro, Block, Model.num_of_qn, "No", Model.p_thread);
   W_Basis[1].qn1 = Model.tot_ele_1;
   W_Basis[1].qn2 = Model.tot_ele_2;
   W_Basis[1].qn3 = Model.tot_sz - 2;
   DMRG_Get_Basis_LLLRRRRL(W_Basis[1], Basis_System, Basis_Enviro, Block, Model.num_of_qn, "No", Model.p_thread);
   
   std::vector<double> SxL(Model.system_size, 0.0)  , SyL(Model.system_size, 0.0);
   std::vector<double> SxC_1(Model.system_size, 0.0), SyC_1(Model.system_size, 0.0);
   std::vector<double> SxC_2(Model.system_size, 0.0), SyC_2(Model.system_size, 0.0);
   
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
   
   DMRG_Transform_Matrix_One(Model.SxC_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxC_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.SxC_1_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SxC_1_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SxC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SxC_1_CF, SxC_1, "SxC_1_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.iSyC_1_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.iSyC_1_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.iSyC_1_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SyC_1_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.iSyC_1_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SyC_1_CF, SyC_1, "SyC_1_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SxC_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxC_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.SxC_2_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SxC_2_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.SxC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SxC_2_CF, SxC_2, "SxC_2_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.iSyC_2_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.iSyC_2_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.iSyC_2_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(SyC_2_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.iSyC_2_On, Mat_Enviro, GS.Vector, Bases_LLLRRRRL, W_Basis, LL_site, RR_site, Model.p_thread);
   Output_Intersite_Values(SyC_2_CF, SyC_2, "SyC_2_CF.txt", LL_site, RR_site, Dmrg_Param, Model);
   
}
