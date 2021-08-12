//
//  Expectation_Values.cpp
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/21.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Expectation_Values(std::vector<double> &GS_Vector, DMRG_LLLRRRRL_Basis &Bases_LLLRRRRL, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Quntum_Number &QN, Model_1D_HUBBARD &Model, DMRG_Param &Dmrg_Param) {
   
   int c1 = (Dmrg_Param.now_sweep == 0);
   int c2 = (Dmrg_Param.tot_sweep - Dmrg_Param.now_sweep >= 2);
   int c3 = (Block_Compo.LL_site + Block_Compo.RR_site + 4 != Model.system_size);
   int c4 = (Block_Compo.LL_site != Block_Compo.RR_site);
   
   if (c1 || c2 || c3 || c4) {
      return;
   }
   
   int LL_site = Block_Compo.LL_site;
   int RR_site = Block_Compo.RR_site;
   
   std::vector<CRS> Mat_System, Mat_Enviro, Mat_CF;
   std::vector<double> Ham, Sz, SxSx, SySy, SzSz, NCUp, NCDown, NC;
   std::vector<double> Sz_CF, NCUp_CF, NCDown_CF, NC_CF, Sx_CF, Sy_CF;
   std::vector<int> Dummy;
   
   DMRG_Transform_Matrix_One(Model.Ham_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.Ham_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(Ham, Mat_System, Model.Ham_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(Ham, "Ham.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.Sz_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.Sz_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(Sz, Mat_System, Model.Sz_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(Sz, "Sz.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.Sz_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(Sz_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.Sz_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(Sz_CF, "Sz_CF.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SxSx_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SxSx_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SxSx, Mat_System, Model.SxSx_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(SxSx, "SxSx.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SySy_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SySy_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SySy, Mat_System, Model.SySy_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(SySy, "SySy.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.SzSz_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.SzSz_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(SzSz, Mat_System, Model.SzSz_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(SzSz, "SzSz.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NCUp_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCUp_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCUp, Mat_System, Model.NCUp_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(NCUp, "NCUp.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NCUp_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCUp_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCUp_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(NCUp_CF, "NCUp_CF.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NCDown_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NCDown_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NCDown, Mat_System, Model.NCDown_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(NCDown, "NCDown.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NCDown_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NCDown_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NCDown_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(NCDown_CF, "NCDown_CF.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.NC_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.NC_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Expectation_Onsite(NC, Mat_System, Model.NC_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Onsite_Values(NC, "NC.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_Two(Model.NC_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(NC_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.NC_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, Block_Compo, Model.p_thread);
   Output_Intersite_Values(NC_CF, "NC_CF.txt", Block_Compo, Dmrg_Param, Model);
   
   std::vector<DMRG_LLLRRRRL_Basis> W_Basis(2);
   DMRG_Quntum_Number W_QN;
   W_QN.qn1_LLLRRRRL = Model.tot_ele;
   W_QN.qn2_LLLRRRRL = Model.tot_sz + 2;
   DMRG_Get_Basis_LLLRRRRL(W_QN, W_Basis[0], Block_Compo, Basis_System, Basis_Enviro, Model.num_of_qn, "No", Model.p_thread);
   W_QN.qn1_LLLRRRRL = Model.tot_ele;
   W_QN.qn2_LLLRRRRL = Model.tot_sz - 2;
   DMRG_Get_Basis_LLLRRRRL(W_QN, W_Basis[1], Block_Compo, Basis_System, Basis_Enviro, Model.num_of_qn, "No", Model.p_thread);
   
   DMRG_Transform_Matrix_One(Model.Sx_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.Sx_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.Sx_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(Sx_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.Sx_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, W_Basis, Block_Compo, Model.p_thread);
   Output_Intersite_Values(Sx_CF, "Sx_CF.txt", Block_Compo, Dmrg_Param, Model);
   
   DMRG_Transform_Matrix_One(Model.iSy_On, Mat_System, Basis_System, LL_site, Model.p_thread);
   DMRG_Transform_Matrix_One(Model.iSy_On, Mat_Enviro, Basis_Enviro, RR_site, Model.p_thread);
   DMRG_Transform_Matrix_Two(Model.iSy_On, Mat_CF, Basis_System, Dummy, "No", LL_site, Model.site_cf_ref, Model.p_thread);
   DMRG_Expectation_Intersite(Sy_CF, Model.site_cf_ref, Mat_CF, Mat_System, Model.iSy_On, Mat_Enviro, GS_Vector, Bases_LLLRRRRL, W_Basis, Block_Compo, Model.p_thread);
   Output_Intersite_Values(Sy_CF, "Sy_CF.txt", Block_Compo, Dmrg_Param, Model);
   
}
