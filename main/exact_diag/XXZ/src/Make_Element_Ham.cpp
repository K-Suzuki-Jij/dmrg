//
//  Make_Element_Ham.cpp
//  1D_XXZ_ED
//
//  Created by Kohei Suzuki on 2020/04/15.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.hpp"

void Make_Element_Ham(long basis, A_Basis_Set &A_Set, Model_1D_XXZ &Model) {
   
   for (int site = 0; site < Model.system_size; site++) {
      A_Set.Local_Basis[site] = ED_Find_Local_Basis(basis, site, Model.dim_onsite);
   }
   
   //Onsite Element
   for (int site = 0; site < Model.system_size; site++) {
      ED_Make_Onsite_Element(basis, site, Model.Ham_On, 1.0, A_Set);
   }
   
   //Intersite Element
   for (int site = 0; site < Model.system_size - 1; site++) {
      ED_Make_Intersite_Element(basis, site, site + 1, Model.Sz_On, Model.Sz_On, Model.J_z     , 1.0, A_Set);
      ED_Make_Intersite_Element(basis, site, site + 1, Model.Sp_On, Model.Sm_On, 0.5*Model.J_xy, 1.0, A_Set);
      ED_Make_Intersite_Element(basis, site, site + 1, Model.Sm_On, Model.Sp_On, 0.5*Model.J_xy, 1.0, A_Set);
   }
   
   if (Model.BC == "PBC") {
      ED_Make_Intersite_Element(basis, 0, Model.system_size - 1, Model.Sz_On, Model.Sz_On, Model.J_z     , 1.0, A_Set);
      ED_Make_Intersite_Element(basis, 0, Model.system_size - 1, Model.Sp_On, Model.Sm_On, 0.5*Model.J_xy, 1.0, A_Set);
      ED_Make_Intersite_Element(basis, 0, Model.system_size - 1, Model.Sm_On, Model.Sp_On, 0.5*Model.J_xy, 1.0, A_Set);
   }
   
   ED_Make_Zero_Element(basis, A_Set);

}
