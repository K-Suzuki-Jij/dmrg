#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "ED.hpp"


void ED_Make_Intersite_Element(long basis, int site1, int site2, CRS &M1, CRS &M2, double coeef, int sign, A_Basis_Set &A_Set) {
   
   if (std::abs(coeef) < A_Set.zero_precision) {
      return;
   }
   
   int local_basis1 = A_Set.Local_Basis[site1];
   int local_basis2 = A_Set.Local_Basis[site2];

   long site_constant1 = A_Set.Site_Constant[site1];
   long site_constant2 = A_Set.Site_Constant[site2];

   for (long i1 = M1.Row[local_basis1]; i1 < M1.Row[local_basis1+1]; i1++) {
      for (long i2 = M2.Row[local_basis2]; i2 < M2.Row[local_basis2+1]; i2++) {
         long a_basis = basis + (M1.Col[i1] - local_basis1)*site_constant1 + (M2.Col[i2] - local_basis2)*site_constant2;
         if (A_Set.Check_Basis.count(a_basis) == 0) {
            A_Set.Check_Basis[a_basis] = (int)A_Set.Basis.size();
            A_Set.Val.push_back(sign*coeef*M1.Val[i1]*M2.Val[i2]);
            A_Set.Basis.push_back(a_basis);
         }
         else {
            A_Set.Val[A_Set.Check_Basis[a_basis]] += sign*coeef*M1.Val[i1]*M2.Val[i2];
         }
         
      }
   }
   
}
