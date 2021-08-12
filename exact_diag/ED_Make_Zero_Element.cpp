#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "ED.hpp"

void ED_Make_Zero_Element(long basis, A_Basis_Set &A_Set) {
   
   if (A_Set.Check_Basis.count(basis) == 0) {
      A_Set.Check_Basis[basis] = (int)A_Set.Basis.size();
      A_Set.Val.push_back(0.0);
      A_Set.Basis.push_back(basis);
   }
   
   
}
