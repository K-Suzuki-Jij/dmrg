//
//  Header.h
//  1D_HUBBARD_DMRG
//
//  Created by Kohei Suzuki on 2020/07/15.
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
#include "Model_1D_HUBBARD.hpp"
#include "DMRG.hpp"
#include "SML.hpp"

struct Block_Operator {
  
   std::vector<CRS> Ham;
   std::vector<CRS> Sz_RE, Sp_RE, Sm_RE;
   std::vector<CRS> CUp_RE, CDown_RE, CUp_D_RE, CDown_D_RE, NC_RE;
   std::vector<CRS> Sz_LE, Sp_LE, Sm_LE;
   std::vector<CRS> CUp_LE, CDown_LE, CUp_D_LE, CDown_D_LE, NC_LE;

};

void DMRG(Model_1D_HUBBARD &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Allocate_Block(Block_Operator &Block, Model_1D_HUBBARD &Model);
void Renormalize(Block_Operator &System, Block_Operator &Enviro, DMRG_Basic_Information &Info, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, Model_1D_HUBBARD &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Get_Operator_Edge(Block_Operator &Block, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_Stored, Model_1D_HUBBARD &Model);
void Find_GS_QN(DMRG_Quntum_Number &QN, DMRG_Block_Composition &Block_Compo, Model_1D_HUBBARD &Model);
void Diagonalize_Ham_LLLRRRRL(Block_Operator &System, Block_Operator &Enviro, DMRG_Basic_Information &Info, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, Model_1D_HUBBARD &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param, DMRG_Time &Time);
void Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, DMRG_LLLRRRRL_Basis &Basis_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Composition &Block_Compo, Model_1D_HUBBARD &Model);
void Make_Elem_Ham_LLLRRRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, std::vector<int> &Ele_LL, std::vector<int> &Ele_On, std::vector<int> &Ele_RR, Model_1D_HUBBARD &Model);
void Renormalize_System(Block_Operator &System, DMRG_LLLR_Basis &Basis_LLLR, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_System, DMRG_T_Mat &T_Mat, Model_1D_HUBBARD &Model);
void Get_Ham_LLLR(Block_Operator &System, DMRG_LLLR_Basis &Basis_LLLR, DMRG_Block_Composition &Block_Compo, std::vector<int> &Ele_LL, Model_1D_HUBBARD &Model, CRS &Ham_LLLR);
void Make_Elem_Ham_LLLR(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLR, Block_Operator &System, std::vector<int> &Ele_LL, Model_1D_HUBBARD &Model);
void Expectation_Values(std::vector<double> &GS_Vector, DMRG_LLLRRRRL_Basis &Bases_LLLRRRRL, DMRG_Block_Composition &Block_Compo, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Quntum_Number &QN, Model_1D_HUBBARD &Model, DMRG_Param &Dmrg_Param);
void Output_Onsite_Values(std::vector<double> &Val, std::string file_name, DMRG_Block_Composition &Block_Compo, DMRG_Param &Dmrg_Param, Model_1D_HUBBARD &Model);
void Output_Intersite_Values(std::vector<double> &Val, std::string file_name, DMRG_Block_Composition &Block_Compo, DMRG_Param &Dmrg_Param, Model_1D_HUBBARD &Model);
#endif /* Header_hpp */
