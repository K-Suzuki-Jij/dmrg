//
//  Header.h
//  1D_AKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/07/08.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#ifndef Header_hpp
#define Header_hpp

#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <algorithm>
#include <ios>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>
#include <cassert>
#include <filesystem>
#include "Model_1D_AKLM.hpp"
#include "DMRG.hpp"
#include "SML.hpp"


struct Block_Operator {
  
   std::vector<CRS> Ham;
   std::vector<CRS> SzL_RE, SpL_RE, SmL_RE;
   std::vector<CRS> SzC_RE, SpC_RE, SmC_RE;
   std::vector<CRS> CUp_RE, CDown_RE, CUp_D_RE, CDown_D_RE, NC_RE;
   std::vector<CRS> SzL_LE, SpL_LE, SmL_LE;
   std::vector<CRS> SzC_LE, SpC_LE, SmC_LE;
   std::vector<CRS> CUp_LE, CDown_LE, CUp_D_LE, CDown_D_LE, NC_LE;

};


void Allocate_Block(Block_Operator &Block, Model_1D_AKLM &Model);
void Allocate_Basis(DMRG_Basis_Stored &Basis, Model_1D_AKLM &Model);
void DMRG(Model_1D_AKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Renormalize(Block_Operator &System, Block_Operator &Enviro, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, DMRG_Basis &Basis, DMRG_Ground_State &GS, Model_1D_AKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Find_GS_QN(int &qn1, int &qn2, DMRG_Block_Information &Block, Model_1D_AKLM &Model);
void Diagonalize_Ham_LLLRRRRL(DMRG_Ground_State &GS, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_AKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param, DMRG_Time &Time);
void Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_AKLM &Model);
void Renormalize_System(Block_Operator &System, DMRG_Basis_LLLR &Basis_LLLR, DMRG_Basis_Stored &Basis_System, int LL_site, DMRG_T_Mat &T_Mat, Model_1D_AKLM &Model);
void Get_Ham_LLLR(Block_Operator &System, DMRG_Basis_LLLR &Basis_LLLR, std::vector<int> &Ele_LL, int LL_site, Model_1D_AKLM &Model, CRS &Ham_LLLR);
void Get_Ham_RRRL(Block_Operator &Enviro, DMRG_Basis_RRRL &Basis_RRRL, std::vector<int> &Ele_RR, int RR_site, Model_1D_AKLM &Model, CRS &Ham_RRRL);
void Get_Ham_LRRL(DMRG_Basis_LRRL &Basis_LRRL, std::vector<int> &Ele_LR, Model_1D_AKLM &Model, CRS &Ham_LRRL);
void Make_Elem_Ham_LLLR(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLR, Block_Operator &System, std::vector<int> &Ele_LL, Model_1D_AKLM &Model);
void Make_Elem_Ham_RRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_RRRL, Block_Operator &Enviro, std::vector<int> &Ele_RR, Model_1D_AKLM &Model);
void Make_Elem_Ham_LRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LRRL, std::vector<int> &Ele_LR, Model_1D_AKLM &Model);
void Expectation_Values(DMRG_Ground_State &GS, DMRG_Basis_LLLRRRRL &Bases_LLLRRRRL, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_AKLM &Model, DMRG_Param &Dmrg_Param);
void Output_Onsite_Values(std::vector<double> &Val, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_AKLM &Model);
void Output_Intersite_Values(std::vector<double> &Val_CF, std::vector<double> &Val_On, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_AKLM &Model);
void Print_Status(DMRG_Param &Dmrg_Param, Model_1D_AKLM &Model);
#endif /* Header_hpp */
