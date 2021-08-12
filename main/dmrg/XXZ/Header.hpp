//
//  Header.hpp
//  1D_XXZ_DMRG
//
//  Created by Kohei Suzuki on 2020/05/18.
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
#include "Model_1D_XXZ.hpp"
#include "DMRG.hpp"
#include "SML.hpp"


struct Block_Operator {
  
   std::vector<CRS> Ham;
   std::vector<CRS> Sz_RE, Sp_RE, Sm_RE;
   std::vector<CRS> Sz_LE, Sp_LE, Sm_LE;

};


void DMRG(Model_1D_XXZ &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Allocate_Block(Block_Operator &Block, Model_1D_XXZ &Model);
void Allocate_Basis(DMRG_Basis_Stored &Basis, Model_1D_XXZ &Model);
void Renormalize(Block_Operator &System, Block_Operator &Enviro, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, DMRG_Basis &Basis, DMRG_Ground_State &GS, Model_1D_XXZ &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Print_Status(DMRG_Param &Dmrg_Param, Model_1D_XXZ &Model);
void Get_Ham_LLLR(Block_Operator &System, DMRG_Basis_LLLR &Basis_LLLR, int LL_site, Model_1D_XXZ &Model, CRS &Ham_LLLR);
void Get_Ham_RRRL(Block_Operator &Enviro, DMRG_Basis_RRRL &Basis_RRRL, int RR_site, Model_1D_XXZ &Model, CRS &Ham_RRRL);
void Make_Elem_Ham_LRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LRRL, Model_1D_XXZ &Model);
void Get_Ham_LRRL(DMRG_Basis_LRRL &Basis_LRRL, Model_1D_XXZ &Model, CRS &Ham_LRRL);
void Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_XXZ &Model);
void Make_Elem_Ham_LLLRRRRL(DMRG_Onsite_Basis &Basis_Onsite, DMRG_A_Basis_Set &A_Basis, DMRG_Basis &Basis, Block_Operator &System, Block_Operator &Enviro, DMRG_Block_Hamiltonian &Ham, Model_1D_XXZ &Model);
void Make_Elem_Ham_LLLR(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LLLR, Block_Operator &System, Model_1D_XXZ &Model);
void Make_Elem_Ham_LRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_LRRL, Model_1D_XXZ &Model);
void Make_Elem_Ham_RRRL(DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, std::vector<int> &Inv_RRRL, Block_Operator &Enviro, Model_1D_XXZ &Model);
void Renormalize_System(Block_Operator &System, DMRG_Basis_LLLR &Basis_LLLR, DMRG_Basis_Stored &Basis_System, int LL_site, DMRG_T_Mat &T_Mat, Model_1D_XXZ &Model);
void Diagonalize_Ham_LLLRRRRL(DMRG_Ground_State &GS, Block_Operator &System, Block_Operator &Enviro, DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_XXZ &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param, DMRG_Time &Time);
void Expectation_Values(DMRG_Ground_State &GS, DMRG_Basis_LLLRRRRL &Bases_LLLRRRRL, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, Model_1D_XXZ &Model, DMRG_Param &Dmrg_Param);
#endif /* Header_hpp */
