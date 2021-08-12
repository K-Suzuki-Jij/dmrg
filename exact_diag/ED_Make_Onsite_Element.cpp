#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "ED.hpp"


void ED_Make_Onsite_Element(long basis, int site, CRS &M, double coeef, A_Basis_Set &A_Set) {
   
   if (std::abs(coeef) < A_Set.zero_precision) {
      return;
   }
   
   int  local_basis   = A_Set.Local_Basis[site];
   long site_constant = A_Set.Site_Constant[site];

   for (long i = M.Row[local_basis]; i < M.Row[local_basis+1]; i++) {
      long a_basis = basis + (M.Col[i] - local_basis)*site_constant;

      if (A_Set.Check_Basis.count(a_basis) == 0) {
         A_Set.Check_Basis[a_basis] = (int)A_Set.Basis.size();
         A_Set.Val.push_back(coeef*M.Val[i]);
         A_Set.Basis.push_back(a_basis);
      }
      else {
         A_Set.Val[A_Set.Check_Basis[a_basis]] += coeef*M.Val[i];
      }
      
   }
   
   
}
