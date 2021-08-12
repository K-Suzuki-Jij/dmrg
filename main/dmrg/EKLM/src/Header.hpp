//
//  Header.hpp
//  1D_EKLM_DMRG
//
//  Created by Kohei Suzuki on 2021/01/01.
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
#include "Model_1D_EKLM.hpp"
#include "DMRG.hpp"
#include "SML.hpp"

struct Block_Operator {
   
   std::vector<CRS> Ham;
   std::vector<std::vector<std::vector<CRS>>> CUp, CDown, CUp_D, CDown_D;
   std::vector<std::vector<std::vector<CRS>>> SpL, SmL, SzL;
   
   std::vector<std::vector<CRS>> NC_Tot;
};

void DMRG(Model_1D_EKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Allocate_Block(Block_Operator &Block, const Model_1D_EKLM &Model);
void Allocate_Basis(DMRG_Basis_Stored &Basis, Model_1D_EKLM &Model);
void Renormalize(Block_Operator &System, Block_Operator &Enviro, DMRG_Basis_Stored &Basis_System, DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, DMRG_Basis &Basis, DMRG_Ground_State &GS, Model_1D_EKLM &Model, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param);
void Find_GS_QN(DMRG_Basis_LLLRRRRL &Basis, const DMRG_Block_Information &Block, const Model_1D_EKLM &Model);
void Make_Elem_Ham(std::string Mat_Type, const DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv, const Block_Operator &System, const Block_Operator &Enviro, const std::vector<int> &Ele, const Model_1D_EKLM &Model);
void Get_Ham_Block(std::string Mat_Type, const Block_Operator &System, const Block_Operator &Enviro, const DMRG_Basis &Basis, const std::vector<int> &Ele, const DMRG_Block_Information Block, const Model_1D_EKLM &Model, CRS &Ham);
void Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, const Block_Operator &System, const Block_Operator &Enviro, const DMRG_Basis &Basis, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, const DMRG_Block_Information &Block, DMRG_Time &Time, const Model_1D_EKLM &Model);
void Renormalize_System(Block_Operator &System, const DMRG_Basis &Basis, DMRG_Basis_Stored &Basis_System, const DMRG_Block_Information &Block, const DMRG_T_Mat &T_Mat, const Model_1D_EKLM &Model);
void Expectation_Values(const DMRG_Ground_State &GS, const DMRG_Basis_LLLRRRRL &Bases_LLLRRRRL, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, const DMRG_Block_Information &Block, const Model_1D_EKLM &Model, const DMRG_Param &Dmrg_Param);
void Output_Onsite_Values(const std::vector<double> &Val, std::string file_name, const DMRG_Block_Information &Block, const DMRG_Param &Dmrg_Param, const Model_1D_EKLM &Model);
#endif /* Header_hpp */
