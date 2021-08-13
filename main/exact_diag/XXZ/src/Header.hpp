//
//  Header.hpp
//  1D_XXZ_ED
//
//  Created by Kohei Suzuki on 2020/04/12.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#ifndef Header_hpp
#define Header_hpp

#include <cmath>
#include <vector>
#include <omp.h>
#include <algorithm>
#include "Model_1D_XXZ.hpp"
#include "SML.hpp"
#include "ED.hpp"

void Get_Bases(Model_1D_XXZ &Model, std::vector<long> &Bases, int target_sz, double &time);
void Get_Ham(Model_1D_XXZ &Model, std::vector<long> &Bases, CRS &Ham, double &time);
void Make_Element_Ham(long basis, A_Basis_Set &A_Set, Model_1D_XXZ &Model);
void Diagonalize(CRS &Ham, std::vector<double> &GS_Vec, Diag_Param &Param, Model_1D_XXZ &Model, double &time);
void Onsite_Expectation_Values(std::vector<double> &Eigen_Vec, std::vector<long> &Bases, Model_1D_XXZ &Model, double &time);
void Intersite_Expectation_Values(std::vector<double> &Eigen_Vec, std::vector<long> &Bases, Model_1D_XXZ &Model, double &time);
#endif /* Header_hpp */

