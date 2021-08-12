//
//  Header.hpp
//  1D_TAKLM_DMRG
//
//  Created by Kohei Suzuki on 2020/08/04.
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
#include <numeric>
#include <filesystem>
#include "Model_1D_TAKLM.hpp"
#include "DMRG.hpp"
#include "SML.hpp"

struct Block_Operator {
  
   std::vector<CRS> Ham;
   std::vector<CRS> SzL_RE, SpL_RE, SmL_RE;
   std::vector<CRS> SzC_1_RE, SpC_1_RE, SmC_1_RE;
   std::vector<CRS> SzC_2_RE, SpC_2_RE, SmC_2_RE;
   std::vector<CRS> CUp_1_RE, CDown_1_RE, CUp_1_D_RE, CDown_1_D_RE, NC_1_RE;
   std::vector<CRS> CUp_2_RE, CDown_2_RE, CUp_2_D_RE, CDown_2_D_RE, NC_2_RE;
   std::vector<CRS> NC_RE;
   std::vector<CRS> SzL_LE, SpL_LE, SmL_LE;
   std::vector<CRS> SzC_1_LE, SpC_1_LE, SmC_1_LE;
   std::vector<CRS> SzC_2_LE, SpC_2_LE, SmC_2_LE;
   std::vector<CRS> CUp_1_LE, CDown_1_LE, CUp_1_D_LE, CDown_1_D_LE, NC_1_LE;
   std::vector<CRS> CUp_2_LE, CDown_2_LE, CUp_2_D_LE, CDown_2_D_LE, NC_2_LE;
   std::vector<CRS> NC_LE;

};

void DMRG(Model_1D_TAKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Allocate_Block(Block_Operator &Block, Model_1D_TAKLM &Model);
void Print_Status(DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model);
void Renormalize(Block_Operator &System, Block_Operator &Enviro, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, DMRG_Basis &Basis, DMRG_Ground_State &GS, Model_1D_TAKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Find_GS_QN(int &qn1, int &qn2, int &qn3, DMRG_Block_Information &Block, Model_1D_TAKLM &Model);
void Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, DMRG_Time &Time, Model_1D_TAKLM &Model);
void Renormalize_System(Block_Operator &System, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, int LL_site, DMRG_T_Mat &T_Mat, Model_1D_TAKLM &Model);
void Get_Ham_Block(std::string Mat_Type, const Block_Operator &Block_Op, const DMRG_Basis &Basis, const std::vector<int> &Ele, int LL_site, int RR_site, const Model_1D_TAKLM &Model, CRS &Ham);
void Make_Elem_Ham(std::string Mat_Type, const DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv, const Block_Operator &Block_Op, const std::vector<int> &Ele, const Model_1D_TAKLM &Model);
void Expectation_Values(DMRG_Ground_State &GS, DMRG_Basis_LLLRRRRL &Bases_LLLRRRRL, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_TAKLM &Model, DMRG_Param &Dmrg_Param);
void Output_Average_Values(std::vector<double> &Val, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model);
void Output_Average_Values(double val, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model);
void Output_Onsite_Values(std::vector<double> &Val, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model);
void Output_Intersite_Values(std::vector<double> &Val_CF, std::vector<double> &Val_On, std::string file_name, int LL_site, int RR_site, DMRG_Param &Dmrg_Param, Model_1D_TAKLM &Model);
void Allocate_Basis(DMRG_Basis_Stored &Basis, Model_1D_TAKLM &Model);
#endif /* Header_hpp */

